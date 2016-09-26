#ifndef LARLITE_MICHELRECODRIVER_CXX
#define LARLITE_MICHELRECODRIVER_CXX

#include "MichelRecoDriver.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/simch.h"
#include "DataFormat/track.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"


namespace larlite {

  MichelRecoDriver::MichelRecoDriver()
    : _hit_tree(nullptr)
  {
    _name="MichelRecoDriver";
    _fout=0;
    _save_clusters=false;
    _avg_gain = 5.3;
    _chchgain.clear();
    // FIX ME: currently set only plane 2 reconstruction
    SetPlane(2);
  }

  bool MichelRecoDriver::initialize() {

    if(_producer.empty()) {
      print(msg::kERROR,__FUNCTION__,"Input data product producer name not specified!");
      throw std::exception();
    }

    _events_info.clear();
    
    _hit_tree = new TTree("_hit_tree","Hit Tree [all hits in event]");
    _hit_tree->Branch("_run"   , &_run    , "run/I");
    _hit_tree->Branch("_subrun", &_subrun , "subrun/I");
    _hit_tree->Branch("_event" , &_event  , "event/I");
    _hit_tree->Branch("_q_v" , "std::vector<double>" , &_q_v);
    _hit_tree->Branch("_w_v" , "std::vector<double>" , &_w_v);
    _hit_tree->Branch("_t_v" , "std::vector<double>" , &_t_v);
    _hit_tree->Branch("_p_v" , "std::vector<double>" , &_p_v);

    _mgr.Initialize();

    _event_time = 0;
    _event_ctr  = 0;

    return true;
  }
  
  bool MichelRecoDriver::analyze(storage_manager* storage) {


    _event_watch.Start();

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    double w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    double t2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    // Get data products
    auto ev_cluster = storage->get_data<event_cluster>(_producer);

    if(!ev_cluster || ev_cluster->empty()) {
      std::cout<<"No cluster found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }

    // get ID information
    _run = storage->get_data<event_cluster>(_producer)->run();
    _subrun = storage->get_data<event_cluster>(_producer)->subrun();
    _event = storage->get_data<event_cluster>(_producer)->event_id();
    
    michel::EventID id;
    id.run    = _run;
    id.subrun = _subrun;
    id.event  = _event;
    
    for (size_t i=0; i < _events_info.size(); i++){
      auto past_event = _events_info[i];
      if ( (past_event.run == id.run) and (past_event.subrun == id.subrun) and (past_event.event == id.event) ){
	std::cout << "We have already scanned this event...skipping..." << std::endl;
	std::cout << "skipping run : " << id.run << "\t subrun : " << id.subrun << "\t event : " << id.event << std::endl;
	return true;
      }
    }

    _events_info.push_back(id);  
    
    if(!ev_cluster || ev_cluster->empty()) return false;

    // Get hits & association
    event_hit* ev_hit = nullptr;
    auto const& hit_ass_set = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());

    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()) {
      std::cout << "No hit found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }
    
    _q_v.clear();
    _w_v.clear();
    _t_v.clear();
    _p_v.clear();

    // Reaching this point means we have something to process. Prepare.
    _mgr.EventReset();
    
    // Tracker for used-hit index
    std::vector< ::michel::HitPt > all_hits_v;
    all_hits_v.reserve(ev_hit->size());

    for(size_t hit_index=0; hit_index<ev_hit->size(); ++hit_index) {
      auto const& h = (*ev_hit)[hit_index];
      // chrage :

      int chan = h.Channel();

      double q = h.Integral();

      double qcalib = q * ( _avg_gain / _chchgain[chan] );

      double w = h.WireID().Wire * w2cm;
      double t = h.PeakTime() * t2cm;

      unsigned int p = h.WireID().Plane;

      if(p == 2){
	_q_v.push_back(q);
	_w_v.push_back(w);
	_t_v.push_back(t);
	_p_v.push_back(p);
      }

      all_hits_v.emplace_back( qcalib , w, t, hit_index, p );
    }
    
    // all hits vector is registered to the manager for future use (search for local hits not clustered)
    _mgr.RegisterAllHits( std::move(all_hits_v) );

    /// set event info
    _mgr.SetEventInfo(id);

    // Loop over clusters & add them to our algorithm manager
    for ( size_t idx=0; idx < hit_ass_set.size(); idx++){

      auto const& hit_ass = hit_ass_set.at(idx);

      // Prepare our hit-list representation    
      std::vector< ::michel::HitPt > michel_cluster;
      michel_cluster.reserve(hit_ass.size());

      // Loop over hits and fill
      for(auto const& hit_index : hit_ass) {

	auto const& h = (*ev_hit)[hit_index];
	
	if(h.WireID().Plane != 2) continue;

	double qcalib = h.Integral() * (_avg_gain / _chchgain[ h.Channel() ] );
	
	michel_cluster.emplace_back( qcalib,
				     h.WireID().Wire * w2cm,
				     h.PeakTime() * t2cm,
				     hit_index,
				     h.WireID().Plane);
      }
      
      // Append a hit-list (cluster) to a manager if not empty
      if(michel_cluster.size() > _minClusSize){
	bool added = _mgr.Append(std::move(michel_cluster),idx);
      }// if above the hit number threshold
    }// for all clusters

    // Now process
    _mgr.Process();

    auto const& michels = _mgr.GetResult();
    if (michels.size() > 0)
      _hit_tree->Fill();
    
    // We may want to save the michels we found as clusters
    // if this option is selected, take the michel hits
    // and group them in "michel" clusters

    if (_save_clusters){

      // create michel cluster object
      auto electron_cluster_v = storage->get_data<event_cluster>("electron");
      auto photon_cluster_v   = storage->get_data<event_cluster>("photon");
      // create association object for clusters to hits
      auto electron_cluster_ass_v = storage->get_data<event_ass>(electron_cluster_v->name());
      auto photon_cluster_ass_v   = storage->get_data<event_ass>(photon_cluster_v->name());
      // create an association for electron cluster to photon clusters
      auto electron_photon_ass_v  = storage->get_data<event_ass>("electron");

      // set event ID through storage manager
      storage->set_id(storage->get_data<event_cluster>(_producer)->run(),
		      storage->get_data<event_cluster>(_producer)->subrun(),
		      storage->get_data<event_cluster>(_producer)->event_id()); 

      // create a vector where to store the cluster -> hit association
      std::vector<std::vector<unsigned int> > electron_clus_hit_ass_v; // michel cluster -> hit
      std::vector<std::vector<unsigned int> > photon_clus_hit_ass_v; // michel cluster -> hit
      // create a vector where to store electron -> photon clusters associations
      std::vector<std::vector<unsigned int> > electron_photon_clus_ass_v; // electron -> photons
      // we need to loop over the michel hits and add them to new clusters

      // get the vector of MichelCluster objects
      auto const& michels = _mgr.GetResult();
      // for each get the hits associated
      for (auto const& michelClus : michels){

	// prepare an empty cluster
	larlite::cluster clus_electron;
	// set start location
	clus_electron.set_start_wire( michelClus._michel._start._w / w2cm, 0. );
	clus_electron.set_start_tick( michelClus._michel._start._t / t2cm, 0. );

	// get the hits (in "michel" notation) for this cluster
	auto const& michel = michelClus._michel; // this is a vector of HitPt
	// michel is a list of HitPt
	// each one has a HitID_t unique ID which we can use
	// to trace it back to the larlite hit it came from
	// the hit index is indeed the position inside the larlite::hit vector
	// we basically need to create an association between the cluster
	// and all the hits in this michel

	// empty vector where to store hits associated for this specific electron cluster
	std::vector<unsigned int> electron_clus_hits;

	// empty vector where to store photon cluster indices associated to the electron cluster
	std::vector<unsigned int> electron_photon_ass;

	// save michel hits

	// 1st electron hits
	double Qelec = 0;
	for (auto const& michel_elec_hit_idx : michel._electron_hit_idx_v){
	  electron_clus_hits.push_back( michel[ michel_elec_hit_idx ]._id );
	  Qelec += michel[ michel_elec_hit_idx ]._q;
	}
	clus_electron.set_summedADC( Qelec, 0, 0);
	electron_clus_hit_ass_v.push_back(electron_clus_hits);
	electron_cluster_v->push_back(clus_electron);

	// 2nd the various photons
	for (auto const& photon_hit_v : michel._photon_clus_v){
	  // empty vector where to store hits associated for this specific photon cluster
	  std::vector<unsigned int> photon_clus_hits;
	  larlite::cluster clus_photon;
	  // charge for this photon
	  double Qphoton =  0;
	  for (auto const& photon_hit_idx : photon_hit_v){
	    photon_clus_hits.push_back( michel[ photon_hit_idx ]._id );
	    Qphoton += michel[ photon_hit_idx ]._q;
	  }
	  clus_photon.set_summedADC( Qelec, 0, 0);
	  photon_clus_hit_ass_v.push_back( photon_clus_hits );
	  photon_cluster_v->push_back(clus_photon);
	  electron_photon_ass.push_back( photon_cluster_v->size() - 1);
	}

	electron_photon_clus_ass_v.push_back( electron_photon_ass );
	
      }// loop over all found michels
      
      // now save the association information
      electron_cluster_ass_v -> set_association( electron_cluster_v->id(), product_id(data::kHit,ev_hit->name()),               electron_clus_hit_ass_v    );
      photon_cluster_ass_v   -> set_association( photon_cluster_v->id(),   product_id(data::kHit,ev_hit->name()),               photon_clus_hit_ass_v      );
      electron_photon_ass_v  -> set_association( electron_cluster_v->id(), product_id(data::kCluster,photon_cluster_v->name()), electron_photon_clus_ass_v );
      
    }// if we should save cluster

    if (_filter_events && (_mgr.GetResult().size() == 0) )
      return false;

    return true;
  }


  bool MichelRecoDriver::finalize() {

    std::cout << "time/event = " << _event_time/_event_ctr * 1.e6 << std::endl;

    _fout->cd();
    _mgr.Finalize(_fout);
    if (_hit_tree)
      _hit_tree->Write();
    return true;
  }


  void MichelRecoDriver::SetChGain(int ch, double g) {

    if (_chchgain.size() >= ch){
      std::cout << "ERROR : channel " << ch << " out of range" << std::endl;
      throw std::exception();
    }
    
    _chchgain[ch] = g;

    return;
  }


}
#endif
