#ifndef LARLITE_MICHELRECODRIVER_CXX
#define LARLITE_MICHELRECODRIVER_CXX

#include "MichelRecoDriver.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/simch.h"

#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/LArProperties.h"


namespace larlite {

  MichelRecoDriver::MichelRecoDriver()
    : _hit_tree(nullptr)
    , _mc_tree(nullptr)
  {
    _name="MichelRecoDriver";
    _fout=0;
    _save_clusters=false;
    _Efield=0.5;
    _use_mc = false;
    // FIX ME: currently set only plane 2 reconstruction
    SetPlane(2);
  }

  bool MichelRecoDriver::initialize() {

    if(_producer.empty()) {
      print(msg::kERROR,__FUNCTION__,"Input data product producer name not specified!");
      throw std::exception();
    }
    
    _hit_tree = new TTree("_hit_tree","Hit Tree [all hits in event]");
    _hit_tree->Branch("_q_v","std::vector<double>",&_q_v);
    _hit_tree->Branch("_w_v","std::vector<double>",&_w_v);
    _hit_tree->Branch("_t_v","std::vector<double>",&_t_v);
    _hit_tree->Branch("_p_v","std::vector<double>",&_p_v);
    _hit_tree->Branch("_run",&_run,"run/I");
    _hit_tree->Branch("_subrun",&_subrun,"subrun/I");
    _hit_tree->Branch("_event",&_event,"event/I");

    _mc_tree = new TTree("_mc_tree","MC comparison TTree");
    _mc_tree->Branch("_run",&_run,"run/I");
    _mc_tree->Branch("_subrun",&_subrun,"subrun/I");
    _mc_tree->Branch("_event",&_event,"event/I");
    
    _mc_tree->Branch("_mc_energy",&_mc_energy,"mc_energy/D");
    _mc_tree->Branch("_reco_energy",&_reco_energy,"reco_energy/D");
    _mc_tree->Branch("_michel_hit_frac","std::vector<double>",&_michel_hit_frac);

    _mgr.Initialize();

    return true;
  }
  
  bool MichelRecoDriver::analyze(storage_manager* storage) {

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    double w2cm = larutil::GeometryUtilities::GetME()->WireToCm();
    // time->cm (accounting for different operating voltages)
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(_Efield,87); // [cm/us]
    // tick width in time
    double tickWidth = 0.5; // [us]
    double t2cm = tickWidth*driftVel;


    // Get data products
    auto ev_cluster = storage->get_data<event_cluster>(_producer);

    // get ID information
    _run = storage->get_data<event_cluster>(_producer)->run();
    _subrun = storage->get_data<event_cluster>(_producer)->subrun();
    _event = storage->get_data<event_cluster>(_producer)->event_id();
    
    michel::EventID id;
    id.run    = _run;
    id.subrun = _subrun;
    id.event  = _event;

    if(!ev_cluster || ev_cluster->empty()) return false;

    // Get hits & association
    event_hit* ev_hit = nullptr;
    auto const& hit_ass_set = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());

    // If ev_hit is null, failed to find data or assocaition (shouldn't happen)
    if(!ev_hit || ev_hit->empty()) return false;

    _q_v.clear();
    _w_v.clear();
    _t_v.clear();
    _p_v.clear();
    
    // Tracker for used-hit index
    std::vector< ::michel::HitPt > all_hits_v;
    all_hits_v.reserve(ev_hit->size());

    for(size_t hit_index=0; hit_index<ev_hit->size(); ++hit_index) {
      auto const& h = (*ev_hit)[hit_index];
      // chrage :
      double q = h.Integral();
      double w = h.WireID().Wire * w2cm;
      double t = (h.PeakTime()-3200) * t2cm;

      unsigned int p = h.WireID().Plane;

      if(p == 2){
	_q_v.push_back(q);
	_w_v.push_back(w);
	_t_v.push_back(t);
      }

      all_hits_v.emplace_back( h.Integral(),
			       w,
			       t,
			       hit_index,
			       p );
    }

    _hit_tree->Fill();
    
    // Reaching this point means we have something to process. Prepare.
    _mgr.EventReset();

    /// set event info
    _mgr.SetEventInfo(id);
        
    // Loop over clusters & add them to our algorithm manager
    for(auto const& hit_ass : hit_ass_set) {

      // Prepare our hit-list representation    
      std::vector< ::michel::HitPt > michel_cluster;
      michel_cluster.reserve(hit_ass.size());

      // Loop over hits and fill
      for(auto const& hit_index : hit_ass) {

	auto const& h = (*ev_hit)[hit_index];
	
	if(h.WireID().Plane != 2) continue;
	
	michel_cluster.emplace_back( h.Integral(),
				     h.WireID().Wire * w2cm,
				     (h.PeakTime() - 3200) * t2cm,
				     hit_index,
				     h.WireID().Plane);
      }
      
      // Append a hit-list (cluster) to a manager if not empty
      if(michel_cluster.size())
	_mgr.Append(std::move(michel_cluster));
    }

    _mgr.RegisterAllHits( std::move(all_hits_v) );
    
    // Now process
    _mgr.Process();

    // We may want to save the michels we found as clusters
    // if this option is selected, take the michel hits
    // and group them in "michel" clusters
    if (_save_clusters){

      auto michel_cluster = storage->get_data<event_cluster>("michel");
      // create association object for clusters
      auto cluster_ass_v = storage->get_data<event_ass>(michel_cluster->name());
      // set event ID through storage manager
      storage->set_id(storage->get_data<event_cluster>(_producer)->run(),
		      storage->get_data<event_cluster>(_producer)->subrun(),
		      storage->get_data<event_cluster>(_producer)->event_id()); 

      // create a vector where to store the cluster -> hit association
      std::vector<std::vector<unsigned int> > clus_hit_ass_v;
      // we need to loop over the michel hits and add them to new clusters

      // get the vector of MichelCluster objects
      auto const& michels = _mgr.GetResult();
      // for each get the hits associated
      for (auto const& michelClus : michels){
	// prepare an empty cluster
	larlite::cluster clus;
	michel_cluster->push_back(clus);
	// get the hits (in "michel" notation) for this cluster
	auto const& michel = michelClus._michel; // this is a vector of HitPt
	// michel is a list of HitPt
	// each one has a HitID_t unique ID which we can use
	// to trace it back to the larlite hit it came from
	// the hit index is indeed the position inside the larlite::hit vector
	// we basically need to create an association between the cluster
	// and all the hits in this michel
	// empty vector where to store hits associated for this specific cluster
	std::vector<unsigned int> clus_hits;
	for (auto const& michel_hit : michel)
	  clus_hits.push_back(michel_hit._id);
	// add this vector to the associations
	if (clus_hits.size() > 3)
	  clus_hit_ass_v.push_back(clus_hits);
	
      }// loop over all found michels
      // now save the association information
      cluster_ass_v->set_association(michel_cluster->id(),product_id(data::kHit,ev_hit->name()),clus_hit_ass_v);
    }// if we should save cluster

    // if we want to compare with mcshowers, we also add feature to backtrack
    if (_use_mc){

      auto ev_mcshower = storage->get_data<event_mcshower>( "mcreco"  );
      auto ev_simch    = storage->get_data<event_simch>   ( "largeant");      

      if(ev_mcshower->size() == 0) { std::cout << "No mcshower found " << "\n"; throw std::exception(); }
      if(ev_simch->size()    == 0) { std::cout << "No simch    found " << "\n"; throw std::exception(); }
 
      bool made_michel = true;

      _reco_energy = 0;
      _mc_energy = -1;

      //********************************
      // Did we make Michel or no, if so get Q, if not move on
      
      auto const& michels = _mgr.GetResult();
      // get cluster's energy
      if (michels.size() != 0) {
	auto const& mich = michels[0]._michel;
	for (auto const& h : mich)
	  _reco_energy += h._q;
      }
      else
	made_michel = false;

      //********************************
      // Is there a Michel or no, if yes i should search simchannel
      // if not don't fill the tree
      
      if(made_michel) {

	//********************************
	// Get the MeV scale energy for this shower

	for(auto const& mcs : *ev_mcshower){
	  if( (mcs.MotherPdgCode() == 13                &&
	       mcs.Process() == "muMinusCaptureAtRest") &&
	      (mcs.DetProfile().E()/mcs.Start().E()  > 0.5
	       || mcs.DetProfile().E() >= 15) )
	    _mc_energy = mcs.DetProfile().E();
	}
	
	//********************************
	// You will probably complain about this but I need to use backtracker
	// how well we can self contain the class in current repo is another story
	
	auto mc_energy_min = 0;  // MeV 
	auto mc_energy_max = 65; // MeV 
      
	std::vector<std::vector<unsigned int> > g4_trackid_v;
	std::vector<unsigned int>               mc_index_v;
	g4_trackid_v.reserve(ev_mcshower->size());
	mc_index_v  .reserve(ev_mcshower->size());
	
	for(size_t mc_index=0; mc_index < ev_mcshower->size(); ++mc_index) {
	  auto const& mcs = (*ev_mcshower)[mc_index];

	  if( (mcs.MotherPdgCode() == 13                &&
	       mcs.Process() == "muMinusCaptureAtRest") &&
	      (mcs.DetProfile().E()/mcs.Start().E()  > 0.5
	       || mcs.DetProfile().E() >= 15) ) { // same criteria as a few lines above, 1/event I hope
	    
	    double energy = mcs.DetProfile().E();
	    
	    std::vector<unsigned int> id_v;
	    id_v.reserve(mcs.DaughterTrackID().size());
	    
	    if( mc_energy_min < energy && energy < mc_energy_max ) {
	      for(auto const& id : mcs.DaughterTrackID()) {
		if(id == mcs.TrackID()) continue;
		id_v.push_back(id);
	      }
	      id_v.push_back(mcs.TrackID());
	      g4_trackid_v.push_back(id_v);
	      mc_index_v.push_back(mc_index);
	    }
	    
	  }
	}

	if(g4_trackid_v.size() == 0) { std::cout << "No tracks found breh \n"; throw std::exception(); } 
	
        try { _BTAlg.BuildMap(g4_trackid_v, *ev_simch, *ev_hit, hit_ass_set); }
	catch(...) { std::cout << "\n ~~..~~ Exception at build map ~~..~~ \n"; }
	
	auto btalgo = _BTAlg.BTAlg();
	
	std::vector<double> hit_frac_michel;
	hit_frac_michel.reserve(ev_hit->size());
	
	for(const auto& h : *ev_hit) {
	  
	  ::btutil::WireRange_t wire_hit(h.Channel(),h.StartTick(),h.EndTick());

	  //offending malloc
	  auto parts = btalgo.MCQ(wire_hit);

	  double michel_part = parts.at(0);
	  double other_part  = parts.at(1);
	  double hit_frac    = michel_part / ( michel_part + other_part );
	  
	  hit_frac_michel.push_back(hit_frac);
	}
	
	//for(const auto& h : hit_frac_michel) { std::cout << " h: " << h << " "; } std::cout << std::endl;	

	if(hit_frac_michel.size() != ev_hit->size()) {
	  std::cout << "I ran backtracker but for some reason, hit_frac_michel didn't come out same size as "
		    << "ev_hit!!!!!\n";
	  throw std::exception();
	}
	
	
	//********************************
	// How do we really know one of our hits is michel, well it has higher fraction of
	// number of electrons from michel than "not"
	std::swap(_michel_hit_frac,hit_frac_michel);
	
      }
      
      _mc_tree->Fill();
      _michel_hit_frac.clear();
      

      
    }
  

  return true;
  }

  bool MichelRecoDriver::finalize() {
    _mgr.Finalize(_fout);
    if (_hit_tree)            _hit_tree->Write();
    if (_mc_tree && _use_mc)  _mc_tree ->Write();
    return true;
  }

}
#endif
