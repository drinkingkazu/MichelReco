#ifndef LARLITE_MICHELRECODRIVER3D_CXX
#define LARLITE_MICHELRECODRIVER3D_CXX

#include "MichelRecoDriver3D.h"

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

  MichelRecoDriver3D::MichelRecoDriver3D()
    : _hit_tree(nullptr)
  {
    _name="MichelRecoDriver3D";
    _fout=0;
    _save_clusters=false;
    _Efield=0.5;
    // FIX ME: currently set only plane 2 reconstruction
    SetPlane(2);
  }

  bool MichelRecoDriver3D::initialize() {

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
  
  bool MichelRecoDriver3D::analyze(storage_manager* storage) {


    _event_watch.Start();

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    double w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    // time->cm (accounting for different operating voltages)
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(_Efield,87); // [cm/us]
    // tick width in time
    double tickWidth = 0.5; // [us]
    double t2cm = tickWidth*driftVel;

    // get ID information
    _run    = storage->get_data<event_cluster>(_producer)->run();
    _subrun = storage->get_data<event_cluster>(_producer)->subrun();
    _event  = storage->get_data<event_cluster>(_producer)->event_id();

    // grab reconstructed tracks
    auto ev_track = storage->get_data<event_track>(_producer);

    if(!ev_track || ev_track->empty()) {
      std::cout<<"No tracks found. Skipping event: "<<storage->event_id()<<std::endl;
      return false;
    }

    // grab reconstructed tracks
    auto out_track = storage->get_data<event_track>("output");

    storage->set_id(_run,_subrun,_event);

    // 1st thing to do is loop through tracks and create GeoAlgo track objects.
    _geoTrj_v.clear();
    for (size_t i=0; i < ev_track->size(); i++){
      // create geotrack from track
      auto const& trk = ev_track->at(i);
      ::geoalgo::Trajectory_t trj;
      if (trk.NumberTrajectoryPoints() < 2) continue;
      for (size_t n=0; n < trk.NumberTrajectoryPoints(); n++){
	auto const& pt = trk.LocationAtPoint(n);
	trj.push_back ( pt );
      }// for all traj points
      _geoTrj_v.push_back(trj);
    }// for all tracks

    std::cout << " added " << _geoTrj_v.size() << " tracks" << std::endl;


    
    // 1st step: check for 3D kink topology in tracks
    std::vector< std::pair<int,int> > trj_pairs_v;
    std::vector<geoalgo::Point_t>   PoCA_v;
    SearchTrajectoryPairs(trj_pairs_v,PoCA_v);
    std::cout << "found " << trj_pairs_v.size() << " good matches" << std::endl;
    for (auto const& trj_pair : trj_pairs_v){
      out_track->emplace_back( ev_track->at( trj_pair.first  ) );
      out_track->emplace_back( ev_track->at( trj_pair.second ) );
    }
			      

    /*
    
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
      double q = h.Integral();

      double w = h.WireID().Wire * w2cm;
      double t = (h.PeakTime()-3200) * t2cm;

      unsigned int p = h.WireID().Plane;

      if(p == 2){
	_q_v.push_back(q);
	_w_v.push_back(w);
	_t_v.push_back(t);
	_p_v.push_back(p);
      }

      all_hits_v.emplace_back( h.Integral(), w, t, hit_index, p );
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
	
	michel_cluster.emplace_back( h.Integral(),
				     h.WireID().Wire * w2cm,
				     (h.PeakTime() - 3200) * t2cm,
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
    
    // loop through michels
    for (auto const& michel : michels){
      if (michel._has_michel){
	// find index with largest number of hits
	size_t max_hits = 0;
	unsigned short largest_clus_idx = 0;
	auto input_clusters = michel.getInputClusterIndex_v();
	std::cout << "reco'd michel from input clusters [ ";
	for (auto idx : input_clusters){
	  std::cout << idx << " ";
	  if (hit_ass_set[idx].size() > max_hits) { 
	    max_hits = hit_ass_set[idx].size();
	    largest_clus_idx = idx;
	  }
	}
	std::cout << "]. Largest cluster index is " << michel.getLargestCluster().first
		  << " with " << michel.getLargestCluster().second << " hits "
		  << " and matches the expected "<< largest_clus_idx << std::endl;

	int startwire = ev_cluster->at(largest_clus_idx).StartWire();
	int endwire   = ev_cluster->at(largest_clus_idx).EndWire();



      }// if there is a michel
    }// for all output MichelClsuters

    // save [run,subrun,event,clus idx] to an output file for saved michels

    // We may want to save the michels we found as clusters
    // if this option is selected, take the michel hits
    // and group them in "michel" clusters

    if (_save_clusters){

      // create michel cluster object
      auto michel_cluster_v = storage->get_data<event_cluster>("michel");
      // create association object for clusters
      auto michel_cluster_ass_v = storage->get_data<event_ass>(michel_cluster_v->name());
      // create muon cluster object
      auto muon_cluster_v = storage->get_data<event_cluster>("muon");
      // create association object for clusters
      auto muon_cluster_ass_v = storage->get_data<event_ass>(muon_cluster_v->name());

      // set event ID through storage manager
      storage->set_id(storage->get_data<event_cluster>(_producer)->run(),
		      storage->get_data<event_cluster>(_producer)->subrun(),
		      storage->get_data<event_cluster>(_producer)->event_id()); 

      // create a vector where to store the cluster -> hit association
      std::vector<std::vector<unsigned int> > michel_clus_hit_ass_v; // michel cluster -> hit
      std::vector<std::vector<unsigned int> > muon_clus_hit_ass_v;   // muon cluster -> hit
      // we need to loop over the michel hits and add them to new clusters

      // get the vector of MichelCluster objects
      auto const& michels = _mgr.GetResult();
      // for each get the hits associated
      for (auto const& michelClus : michels){

	// prepare an empty cluster
	larlite::cluster clus_michel;
	michel_cluster_v->push_back(clus_michel);
	// get the hits (in "michel" notation) for this cluster
	auto const& michel = michelClus._michel; // this is a vector of HitPt
	// michel is a list of HitPt
	// each one has a HitID_t unique ID which we can use
	// to trace it back to the larlite hit it came from
	// the hit index is indeed the position inside the larlite::hit vector
	// we basically need to create an association between the cluster
	// and all the hits in this michel

	// empty vector where to store hits associated for this specific michel cluster
	std::vector<unsigned int> michel_clus_hits;

	// additionally get hit range in time/wire
	int wire_min = 8256;
	int wire_max = 0;
	double tick_min = 9600;
	double tick_max = 0;

	for (auto const& michel_hit : michel){
	  michel_clus_hits.push_back(michel_hit._id);
	  auto const& hit = ev_hit->at(michel_hit._id);
	  auto const& hit_tick = hit.PeakTime();
	  if (hit_tick > tick_max) { tick_max = hit_tick; }
	  if (hit_tick < tick_min) { tick_min = hit_tick; }
	  auto const& hit_wire = hit.WireID().Wire;
	  if (hit_wire > wire_max) { wire_max = hit_wire; }
	  if (hit_wire < wire_min) { wire_min = hit_wire; }
	}

	std::cout << "Reconstructed Michel w/ " << michel.size() << " hits." << std::endl
		  << "w/ tick-range [" << tick_min << ", " << tick_max << "]" << std::endl
		  << "w/ wire-range [" << wire_min << ", " << wire_max << "]" << std::endl
		  << std::endl;

	// add this vector to the associations
	if (michel_clus_hits.size() > 3){
	  michel_clus_hit_ass_v.push_back(michel_clus_hits);
	  // if we want to save the michel hits,
	  // do the same for the muon hits
	  // prepare an empty cluster
	  larlite::cluster clus_muon;
	  muon_cluster_v->push_back(clus_muon);
	  auto const& muon = michelClus._hits;
	  // empty vector where to store hits associated for this specific muon cluster
	  std::vector<unsigned int> muon_clus_hits;
	  // get the boundary position
	  auto const& boundary = michelClus._boundary;
	  // is the michel forwards or backwards?
	  auto const& forwards = michelClus._forward;
	  if (forwards){			
	    for (size_t n=0; n < boundary; n++)
	      muon_clus_hits.push_back(muon[n]._id);
	  }
	  else{
	    for (size_t n=boundary; n < muon.size(); n++)
	      muon_clus_hits.push_back(muon[n]._id);
	  } 
	  // add the association information
	  muon_clus_hit_ass_v.push_back(muon_clus_hits);
	  
	}// if there are at least 3 hits in the michel cluster
	
      }// loop over all found michels
      // now save the association information
      michel_cluster_ass_v->set_association(michel_cluster_v->id(),product_id(data::kHit,ev_hit->name()),michel_clus_hit_ass_v);
      muon_cluster_ass_v->set_association(muon_cluster_v->id(),product_id(data::kHit,ev_hit->name()),muon_clus_hit_ass_v);
      
    }// if we should save cluster

    if (_filter_events && (_mgr.GetResult().size() == 0) )
      return false;

    */

    return true;
  }


  bool MichelRecoDriver3D::finalize() {

    std::cout << "time/event = " << _event_time/_event_ctr * 1.e6 << std::endl;


    auto ts = ::larutil::TimeService::GetME();
    std::cout << "TPCTick -> TDC: 1000 -> " << ts->TPCTick2TDC(1000) << std::endl;
    std::cout << "TPCTick -> TDC: 2000 -> " << ts->TPCTick2TDC(2000) << std::endl;

    _fout->cd();
    _mgr.Finalize(_fout);
    if (_hit_tree)
      _hit_tree->Write();
    return true;
  }


  void MichelRecoDriver3D::SearchTrajectoryPairs(std::vector< std::pair<int,int> >& trj_pairs_v,
						 std::vector<geoalgo::Point_t>& PoCA_v) {
    
    trj_pairs_v.clear();
    PoCA_v.clear();

    for (size_t n=0; n < _geoTrj_v.size(); n++){

      auto const& trj1 = _geoTrj_v[n];

      if (trj1.Length() < 4) continue;

      for (size_t m=n+1; m < _geoTrj_v.size(); m++){

	auto const& trj2 = _geoTrj_v[m];

	if (trj2.Length() < 4) continue;

	geoalgo::Point_t pt1;
	geoalgo::Point_t pt2;

	double sqDistMin = 1000.;
	double sqDist1 = trj1.at(0).Dist( trj2.at(0) );
	if (sqDist1 < sqDistMin) { sqDistMin = sqDist1; pt1 = trj1.at(0); pt2 = trj2.at(0); }
	double sqDist2 = trj1.at( trj1.size()-1 ).Dist( trj2.at(0) );
	if (sqDist2 < sqDistMin) { sqDistMin = sqDist2; pt1 = trj1.at( trj1.size()-1 ); pt2 = trj2.at(0); }
	double sqDist3 = trj1.at(0).Dist( trj2.at( trj2.size()-1 ) );
	if (sqDist1 < sqDistMin) { sqDistMin = sqDist3; pt1 = trj1.at(0); pt2 = trj2.at( trj2.size()-1 ); }
	double sqDist4 = trj1.at( trj1.size()-1 ).Dist( trj2.at( trj2.size()-1 ) );
	if (sqDist1 < sqDistMin) { sqDistMin = sqDist4; pt1 = trj1.at( trj1.size()-1 ); pt2 = trj2.at( trj2.size()-1 ); }

	if (sqDistMin < 5)
	  trj_pairs_v.push_back( std::make_pair(n,m) );

      }// for all trajectories, loop 2
    }// for all trajectories, loop 1

    return;

  }

}
#endif