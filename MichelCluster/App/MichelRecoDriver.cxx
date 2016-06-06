#ifndef LARLITE_MICHELRECODRIVER_CXX
#define LARLITE_MICHELRECODRIVER_CXX

#include "MichelRecoDriver.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/simch.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/track.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"


namespace larlite {

  MichelRecoDriver::MichelRecoDriver()
    : _hit_tree(nullptr)
    , _MIP_tree(nullptr)
    , _mc_tree(nullptr)
    , _mc_hit_tree(nullptr)
  {
    _name="MichelRecoDriver";
    _fout=0;
    _save_clusters=false;
    _Efield=0.5;
    _use_mc = false;
    _out_txt_file.open("input_cluster_info_file.txt");
    _debug_mcq = false;
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



    _mc_hit_tree = new TTree("_mc_hit_tree","MC Hit Information");
    _mc_hit_tree->Branch("_hit_integral",&_hit_integral,"hit_integral/D");
    _mc_hit_tree->Branch("_hit_mc_q",&_hit_mc_q,"hit_mc_q/D");

    // MC Michel & Muon Information (& Simchannel stuff)
    _mc_tree = new TTree("_mc_tree","MC comparison TTree");
    // event information
    _mc_tree->Branch("_run"   , &_run    , "run/I");
    _mc_tree->Branch("_subrun", &_subrun , "subrun/I");
    _mc_tree->Branch("_event" , &_event  , "event/I");
    // michel MC information
    _mc_tree->Branch("_mc_x",  &_mc_x,  "mc_x/D" );
    _mc_tree->Branch("_mc_y",  &_mc_y,  "mc_y/D" );
    _mc_tree->Branch("_mc_z",  &_mc_z,  "mc_z/D" );
    _mc_tree->Branch("_mc_px", &_mc_px, "mc_px/D");
    _mc_tree->Branch("_mc_py", &_mc_py, "mc_py/D");
    _mc_tree->Branch("_mc_pz", &_mc_pz, "mc_pz/D");
    _mc_tree->Branch("_mc_yplane_angle",&_mc_yplane_angle,"mc_yplane_angle/D");
    // muon MC information
    _mc_tree->Branch("_mu_x",  &_mu_x,  "mu_x/D" );
    _mc_tree->Branch("_mu_y",  &_mu_y,  "mu_y/D" );
    _mc_tree->Branch("_mu_z",  &_mu_z,  "mu_z/D" );
    _mc_tree->Branch("_mu_px", &_mu_px, "mu_px/D");
    _mc_tree->Branch("_mu_py", &_mu_py, "mu_py/D");
    _mc_tree->Branch("_mu_pz", &_mu_pz, "mu_pz/D");
    _mc_tree->Branch("_mu_yplane_angle",&_mu_yplane_angle,"mu_yplane_angle/D");
    // 3D dot product between michel and muon direction
    _mc_tree->Branch("_mu_michel_3dangle",&_mu_michel_3dangle,"mu_michel_3dangle/D");
    // Reconstructed energy & Simch stuff
    _mc_tree->Branch("_mc_energy"      , &_mc_energy           , "mc_energy/D");
    _mc_tree->Branch("_reco_energy"    , &_reco_energy         , "reco_energy/D");
    _mc_tree->Branch("_michel_hit_fracReco", "std::vector<double>" , &_michel_hit_fracReco);
    _mc_tree->Branch("_michel_hit_QtotReco", "std::vector<double>" , &_michel_hit_QtotReco);
    _mc_tree->Branch("_michel_hit_idxReco",  "std::vector<double>" , &_michel_hit_idxReco );
    _mc_tree->Branch("_michel_hit_fracMC", "std::vector<double>" , &_michel_hit_fracMC);
    _mc_tree->Branch("_michel_hit_QtotMC", "std::vector<double>" , &_michel_hit_QtotMC);
    _mc_tree->Branch("_QMichelMC", &_QMichelMC, "QMichelMC/D");
    _mc_tree->Branch("_QMichelRecoSimch_all", &_QMichelRecoSimch_all, "QMichelRecoSimch_all/D");
    _mc_tree->Branch("_QMichelRecoSimch_shr", &_QMichelRecoSimch_shr, "QMichelRecoSimch_shr/D");
    _mc_tree->Branch("_QMichelShowerMCSimch_all", &_QMichelShowerMCSimch_all, "QMichelShowerMCSimch_all/D");
    _mc_tree->Branch("_QMichelShowerMCSimch_shr", &_QMichelShowerMCSimch_shr, "QMichelShowerMCSimch_shr/D");
    _mc_tree->Branch("_QMichelPartMCSimch_all", &_QMichelPartMCSimch_all, "QMichelPartMCSimch_all/D");
    _mc_tree->Branch("_QMichelPartMCSimch_shr", &_QMichelPartMCSimch_shr, "QMichelPartMCSimch_shr/D");
    _mc_tree->Branch("_f_RecoHitsQ_fromMichelSimch",&_f_RecoHitsQ_fromMichelSimch,"f_RecoHitsQ_fromMichelSimch/D");
    _mgr.Initialize();

    // write a tree that stores the hit charge for every hit in the event
    _MIP_tree = new TTree("_MIP_tree","MIP Tree");
    _MIP_tree->Branch("_hit_charge",&_hit_charge,"hit_charge/D");
    

    _event_time = 0;
    _event_ctr  = 0;

    return true;
  }
  
  bool MichelRecoDriver::analyze(storage_manager* storage) {


    _event_watch.Start();

    _QMichelRecoSimch_all     = 0.;
    _QMichelRecoSimch_shr     = 0.;
    _QMichelPartMCSimch_all   = 0.;
    _QMichelPartMCSimch_shr   = 0.;
    _QMichelShowerMCSimch_all = 0.;
    _QMichelShowerMCSimch_shr = 0.;
    _michel_hit_fracReco.clear();
    _michel_hit_QtotReco.clear();
    _michel_hit_idxReco.clear();
    _michel_hit_fracMC.clear();
    _michel_hit_QtotMC.clear();


    _wiremap.clear();

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    double w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    // time->cm (accounting for different operating voltages)
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(_Efield,87); // [cm/us]
    // tick width in time
    double tickWidth = 0.5; // [us]
    double t2cm = tickWidth*driftVel;

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

    // get the tracks associated to the hits
    event_track* ev_track = nullptr;
    auto const& track_ass_set = storage->find_one_ass(ev_hit->id(), ev_track, "pandoraCosmicKHit");
    
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

	_out_txt_file << _run << " " << _subrun << " " << _event << " " << largest_clus_idx << " " << startwire << " " << endwire << "\n";

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
      // create michel particle object
      auto michel_particle_v = storage->get_data<event_pfpart>("michel");
      // create association object for pfpart
      auto michel_particle_michel_clus_ass_v = storage->get_data<event_ass>("michelkClustermichel");
      auto michel_particle_input_clus_ass_v = storage->get_data<event_ass>("michelkClusterinput");
      auto michel_particle_input_track_ass_v = storage->get_data<event_ass>("michelkTrackinput");

      // set event ID through storage manager
      storage->set_id(storage->get_data<event_cluster>(_producer)->run(),
		      storage->get_data<event_cluster>(_producer)->subrun(),
		      storage->get_data<event_cluster>(_producer)->event_id()); 

      // create a vector where to store the cluster -> hit association
      std::vector<std::vector<unsigned int> > michel_clus_hit_ass_v; // michel cluster -> hit
      std::vector<std::vector<unsigned int> > muon_clus_hit_ass_v;   // muon cluster -> hit
      std::vector<std::vector<unsigned int> > michel_part_michel_clus_ass_v; // michel pfpart -> michel cluster
      std::vector<std::vector<unsigned int> > michel_part_input_clus_ass_v;  // michel pfpart -> input cluster
      std::vector<std::vector<unsigned int> > michel_part_input_track_ass_v; // michel pfpart -> input track
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
	for (auto const& michel_hit : michel)
	  michel_clus_hits.push_back(michel_hit._id);
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

	  // create a PFParticle for the michel
	  larlite::pfpart michel_part(11,0,0,std::vector<size_t>());
	  michel_particle_v->push_back(michel_part);

	  // find the 3dtrack associated to the input cluster
	  // 1) hits associated to input cluster
	  auto hits_associated_to_input_cluster = hit_ass_set[ michelClus.getLargestCluster().first ];
	  // 2) track associated to these hits
	  unsigned int trk_idx = 0;
	  bool found = false;
	  for (auto const& hit_idx : hits_associated_to_input_cluster){
	    for (size_t h_idx = 0; h_idx < track_ass_set.size(); h_idx++){
	      if (hit_idx == h_idx){
		if (track_ass_set[h_idx].size() == 0) continue;
		trk_idx = track_ass_set[h_idx][0];
		found = true;
		break;
	      }
	    }
	    if (found) { break; }
	  }
	    
	  
	  // keep track of the association to the input cluster
	  michel_part_input_clus_ass_v.push_back( std::vector<unsigned int>(1, michelClus.getLargestCluster().first) );
	  // keep track of the association to the michel cluster
	  michel_part_michel_clus_ass_v.push_back( std::vector<unsigned int>(1, ((unsigned int)michel_cluster_v->size()) - 1) );
	  // keep track of the association to the input track
	  if (found){
	    michel_part_input_track_ass_v.push_back( std::vector<unsigned int>(1, trk_idx ) );
	  }

	}// if there are at least 3 hits in the michel cluster
	
      }// loop over all found michels
      // now save the association information
      michel_cluster_ass_v->set_association(michel_cluster_v->id(),product_id(data::kHit,ev_hit->name()),michel_clus_hit_ass_v);
      muon_cluster_ass_v->set_association(muon_cluster_v->id(),product_id(data::kHit,ev_hit->name()),muon_clus_hit_ass_v);
      michel_particle_michel_clus_ass_v->set_association( michel_particle_v->id(),
							  product_id(data::kCluster,michel_cluster_v->name()),
							  michel_part_michel_clus_ass_v);
      michel_particle_input_clus_ass_v->set_association ( michel_particle_v->id(),
							  product_id(data::kCluster,ev_cluster->name()),
							  michel_part_input_clus_ass_v);
      //std::cout << "ev track name : " << ev_track->name() << std::endl;
      //std::cout << "ass vector size : " << michel_part_input_track_ass_v.size() << std::endl;
      //std::cout << "event : " << _event << "\t subrun : " << _subrun << "\t run : " << _run << std::endl;
      michel_particle_input_track_ass_v->set_association( michel_particle_v->id(),
							  product_id(data::kTrack,ev_track->name()),
							  michel_part_input_track_ass_v);
    }// if we should save cluster

    // if we want to compare with mcshowers, we also add feature to backtrack
    if (_use_mc){



      auto ev_mcshower = storage->get_data<event_mcshower>( "mcreco"   );
      auto ev_mctrack  = storage->get_data<event_mctrack> ( "mcreco"   );
      auto ev_simch    = storage->get_data<event_simch>   ( "largeant" );

      //if(ev_simch->size() == 0) { std::cout << "No Simchannel found " << "\n"; throw std::exception(); }
 
      bool shower_exists = false;

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
      
      //********************************
      // Get the MeV scale energy for this shower
      
      for(auto const& mcs : *ev_mcshower){
	
	if ( ( ( mcs.MotherPdgCode() == 13 ) || ( mcs.MotherPdgCode() == -13 ) ) &&              // parent is a muon
	     ( ( mcs.Process() == "muMinusCaptureAtRest" ) || ( mcs.Process() == "Decay" ) ) &&  // process is consistent with michel decay
	     ( mcs.DetProfile().E()/mcs.Start().E()  > 0.5 ) ) {                                 // at least 50% of energy is contained
	  
	  
	  shower_exists = true;
	}
      }
      
      //********************************
      // If there is an MC shower in this event, backtrack reconstructed
      // hits and find charge deposition

      bool reco_michel_exists = true;
      if (_mgr.GetResult().size() == 0)
	reco_michel_exists = false;
      if (reco_michel_exists){
	if (_mgr.GetResult()[0]._michel.size() == 0)
	  reco_michel_exists = false;
      }
	  
      
      if(shower_exists){// && reco_michel_exists ) { 
	
	//********************************
	// You will probably complain about this but I need to use backtracker
	// how well we can self contain the class in current repo is another story
	
	auto mc_energy_min = 0;  // MeV 
	auto mc_energy_max = 65; // MeV 

	// gather up track IDs for all particles in the Michel EM Shower
	std::vector<std::vector<unsigned int> > g4_trackid_v_shr;
	// use only the Michel particle, ignoring radiated photons
	std::vector<std::vector<unsigned int> > g4_trackid_v_part;
	
	for(size_t mc_index=0; mc_index < ev_mcshower->size(); ++mc_index) {
	  auto const& mcs = (*ev_mcshower)[mc_index];
	  
	  
	  if ( ( ( mcs.MotherPdgCode() == 13 ) || ( mcs.MotherPdgCode() == -13 ) ) &&              // parent is a muon
	       ( ( mcs.Process() == "muMinusCaptureAtRest" ) || ( mcs.Process() == "Decay" ) ) &&  // process is consistent with michel decay
	       ( mcs.DetProfile().E()/mcs.Start().E()  > 0.5 ) ) {                                 // at least 50% of energy is contained
	    
	    _mc_energy = mcs.DetProfile().E();
	    _QMichelMC = mcs.Charge(2);
	    // calculate lifetime correction
	    _mc_x  = mcs.Start().X();
	    _mc_y  = mcs.Start().Y();
	    _mc_z  = mcs.Start().Z();
	    _mc_px = mcs.Start().Px();
	    _mc_py = mcs.Start().Py();
	    _mc_pz = mcs.Start().Pz();
	    double mag = sqrt( (_mc_px*_mc_px) + (_mc_py*_mc_py) + (_mc_pz*_mc_pz) );	    
	    _mc_px /= mag;
	    _mc_py /= mag;
	    _mc_pz /= mag;
	    // 2D angle on Y plane view (in x-z space)
	    _mc_yplane_angle = _mc_pz / sqrt ( _mc_px*_mc_px + _mc_pz*_mc_pz );

	    // find the MCTrack that produced this Michel electron
	    for (auto const& mct : *ev_mctrack){
	      if ( (mcs.MotherTrackID() == mct.TrackID()) and ( abs(mct.PdgCode()) == 13 ) ){
		// if only 1 step: continue
		if (mct.size() < 2)
		  continue;
		_mu_energy  = mct[0].E();
		_mu_px = mct[0].Px();
		_mu_py = mct[0].Py();
		_mu_pz = mct[0].Pz();
		mag = sqrt( (_mu_px*_mu_px) + (_mu_py*_mu_py) + (_mu_pz*_mu_pz) );
		_mu_px /= mag;
		_mu_py /= mag;
		_mu_pz /= mag;
		// 2D angle on Y plane view (in x-z space)
		_mu_yplane_angle = _mu_pz / sqrt ( _mu_px*_mu_px + _mu_pz*_mu_pz );
		// muon - michel 3D angle
		_mu_michel_3dangle = _mu_px*_mc_px + _mu_py*_mc_py + _mu_pz*_mc_pz;
	      }// if we found the parent muon
	    }// for all mctracks
	    
	    double energy = mcs.DetProfile().E();
	    
	    std::vector<unsigned int> id_v;
	    id_v.reserve(mcs.DaughterTrackID().size());
	    
	    if( mc_energy_min < energy && energy < mc_energy_max ) {
	      for(auto const& id : mcs.DaughterTrackID()) {
		if(id == mcs.TrackID()) continue;
		id_v.push_back(id);
	      }
	      id_v.push_back(mcs.TrackID());
	      g4_trackid_v_shr.push_back(id_v);
	      // only for michel particle and no photons
	      std::vector<unsigned int> id_v_part = {mcs.TrackID()};
	      g4_trackid_v_part.push_back(id_v_part);
	    }
	  }
	}
	
	if(g4_trackid_v_shr.size() == 0)
	  std::cout << "No track IDs found for backtracker !!!" << std::endl;
	
	else{
	  
	  try { _BTAlgShower.BuildMap(g4_trackid_v_shr, *ev_simch, *ev_hit, hit_ass_set); }
	  catch(...) { std::cout << "\n Exception at building BTAlg map...\n"; }
	  try { _BTAlgPart.BuildMap(g4_trackid_v_part, *ev_simch, *ev_hit, hit_ass_set); }
	  catch(...) { std::cout << "\n Exception at building BTAlg map...\n"; }
	  
	  auto btalgoShower = _BTAlgShower.BTAlg();
	  auto btalgoPart   = _BTAlgPart.BTAlg();
	  
	  std::vector<double> hit_frac_michel;
	  std::vector<double> hit_Qtot_michel;
	  
	  // hits used for michel
	  int michel_hits = 0;
	  
	  for (size_t hit_index = 0; hit_index < ev_hit->size(); hit_index++){
	    
	    const auto& h = ev_hit->at(hit_index);
	    
	    if(h.WireID().Plane != 2)
	      continue;
	    
	    //::btutil::WireRange_t wire_hit(h.Channel(),h.StartTick()+3050,h.EndTick()+3050);
	    ::btutil::WireRange_t wire_hit(h.Channel(),h.PeakTime()-3*h.RMS()+3050,h.PeakTime()+3*h.RMS()+3050);
	    
	    // check if this wire-range has already been used
	    auto wire_ranges = getUnUsedWireRange(wire_hit);
	    
	    // loop through all output wire-ranges
	    for (auto const& wire_range : wire_ranges){
	      
	      //if the wire-ranges are empty, continue
	      if (wire_range.start == wire_range.end)
		continue;
	      
	      if (_debug_mcq)
		std::cout << "hit wire : " << h.Channel() << " w/ range [" << wire_range.start << ", " << wire_range.end << "]" << " with integral : " << h.Integral() << std::endl;
	      
	      // offending malloc if _BTAlg is not class member...
	      auto partsShower = btalgoShower.MCQ(wire_range);
	      auto partsPart   = btalgoPart.MCQ(wire_range);
	      
	      // calculate hit MC information
	      _hit_integral = h.Integral();
	      _hit_mc_q     = partsShower.at(0)+partsShower.at(1);
	      _mc_hit_tree->Fill();
	      
	      
	      double michel_part = partsShower.at(0);
	      double other_part  = partsShower.at(1);
	      double hit_frac    = michel_part / ( michel_part + other_part );
	      if (_debug_mcq)
		std::cout << "hit " << hit_index << " has " << michel_part+other_part << " and " << hit_frac << " michel contribution..." << std::endl;
	      if (hit_frac > 0.1){
		_QMichelShowerMCSimch_shr += michel_part;
		_QMichelShowerMCSimch_all += (michel_part+other_part);
		_michel_hit_QtotMC.push_back(michel_part);
		_michel_hit_fracMC.push_back(hit_frac);
	      }
	      
	      // if no michel is reconstructed, skip this part
	      if (reco_michel_exists){
		// check if this hit has been added to the michel cluster...
		if (_mgr.GetResult().size() == 0) continue;
		auto michel = _mgr.GetResult()[0]._michel;
		for (auto const& h : michel){
		  if (h._id == hit_index){
		    if (_debug_mcq)
		      std::cout << "...this hit is michel" << std::endl;
		    michel_hits += 1;
		    _QMichelRecoSimch_all += (michel_part+other_part);
		    _QMichelRecoSimch_shr += michel_part;
		    _michel_hit_QtotReco.push_back(michel_part);
		    _michel_hit_fracReco.push_back(hit_frac);
		    _michel_hit_idxReco.push_back(h._id);
		  }
		}// for all michel hits
	      }// if a reconstructed michel exists
	      
	      // now add charge from michel original electron only
	      michel_part = partsPart.at(0);
	      other_part  = partsPart.at(1);
	      hit_frac    = michel_part / ( michel_part + other_part );
	      if (hit_frac > 0.1){
		_QMichelPartMCSimch_shr += michel_part;
		_QMichelPartMCSimch_all += michel_part + other_part;
	      }
	      
	    } // for all wire-ranges
	  } // for all hits

	  _f_RecoHitsQ_fromMichelSimch = _QMichelRecoSimch_shr / _QMichelRecoSimch_all;
	  
	  if (_debug_mcq){
	    std::cout << std::endl
		      << "QMichel Reconstructed (all) (simch) : " << _QMichelRecoSimch_all << std::endl
		      << "QMichel Reconstructed (shr) (simch) : " << _QMichelRecoSimch_shr << std::endl
		      << "QMichel Shower MC (all) (simch)     : " << _QMichelShowerMCSimch_all << std::endl
		      << "QMichel Shower MC (shr) (simch)     : " << _QMichelShowerMCSimch_shr << std::endl
		      << "QMichel Part MC (all) (simch)       : " << _QMichelPartMCSimch_all << std::endl
		      << "QMichel Part MC (shr) (simch)       : " << _QMichelPartMCSimch_shr << std::endl
		      << "Frac of Q actually from the michel  : " << _f_RecoHitsQ_fromMichelSimch << std::endl
		      << "Frac of  MC Shower reco'd in hits   : " << _QMichelRecoSimch_all/_QMichelShowerMCSimch_shr << std::endl
		      << "Frac of  MC Part   reco'd in hits   : " << _QMichelRecoSimch_all/_QMichelPartMCSimch_shr << std::endl
		      << std::endl << std::endl;
	  }
	  
	}// if trackIDs are found for the back-tracker
	//********************************
	// How do we really know one of our hits is michel, well it has higher fraction of
	// number of electrons from michel than "not", go ahead and swap for writeout, if this vector
	// is empty in TTree then there was no MC shower
	
	
      } // if we found an mcshower

      _mc_tree->Fill();
    }

    _event_time += _event_watch.RealTime();    
    _event_ctr  += 1;  

  return true;
  }

  bool MichelRecoDriver::finalize() {

    _out_txt_file.close();

    std::cout << "time/event = " << _event_time/_event_ctr * 1.e6 << std::endl;


    auto ts = ::larutil::TimeService::GetME();
    std::cout << "TPCTick -> TDC: 1000 -> " << ts->TPCTick2TDC(1000) << std::endl;
    std::cout << "TPCTick -> TDC: 2000 -> " << ts->TPCTick2TDC(2000) << std::endl;

    _fout->cd();
    _mgr.Finalize(_fout);
    if (_hit_tree)
      _hit_tree->Write();
    if (_mc_tree && _use_mc)
      _mc_tree ->Write();
    if (_mc_hit_tree)
      _mc_hit_tree->Write();
    //if (_MIP_tree)
    //  _MIP_tree->Write();
    return true;
  }


  std::vector<btutil::WireRange_t> MichelRecoDriver::getUnUsedWireRange(const btutil::WireRange_t& wirerange){

    //std::cout << "searching for already used wire-range @ chan " << wirerange.ch << " in range [" << wirerange.start << ", " << wirerange.end << "]" << std::endl;

    // return a vector of WireRanges
    std::vector<btutil::WireRange_t> out_ranges;

    // is this wire in the map?
    // YES
    if ( _wiremap.find( wirerange.ch ) == _wiremap.end() ){
      // no -> add to the map but don't edit anything
      //std::cout << "\t this channel has not been added..." << std::endl;
      std::pair<double,double> tickrange = std::make_pair(wirerange.start,wirerange.end);
      std::vector< std::pair<double,double> > tickranges = {tickrange};
      _wiremap[ wirerange.ch ] = tickranges;
      // return the same exact wire-range
      //std::cout << "\t adding wire-range @ chan "<< wirerange.ch << " in range [" << wirerange.start << ", " << wirerange.end << "]" << std::endl;
      out_ranges.push_back( wirerange );
    }
    // NO
    else{
      // the channel is there, but is the same time-range?
      auto& wire_ranges = _wiremap[ wirerange.ch ];
      // find if these ticks are included
      // has a wire-range overlap been found?
      bool found = false;
      for (size_t i=0; i < wire_ranges.size(); i++){
	auto wire_range = wire_ranges[i];
	// if wirerange is before 
	if ( wirerange.end < wire_range.first )
	  continue;
	// if wirarange is after
	if ( wirerange.start > wire_range.second )
	  continue;
	// if we made it this far -> there is overlap
	// find the overlap region
	// the two completely overlap, return an empty wire-range
	if ( (wirerange.start >= wire_range.first) and (wirerange.end <= wire_range.second) ){
	  found = true;
	  btutil::WireRange_t newrange( wirerange.ch, 0, 0);
	  //std::cout << "\t adding wire-range @ chan "<< newrange.ch << " in range [" << newrange.start << ", " << newrange.end << "]" << std::endl;
	  out_ranges.push_back ( newrange );
	}
	// if we overflow on the low-end
	if ( wirerange.start < wire_range.first ){
	  found = true;
	  btutil::WireRange_t newrange( wirerange.ch, wirerange.start, wire_range.first);
	  //std::cout << "\t adding wire-range @ chan "<< newrange.ch << " in range [" << newrange.start << ", " << newrange.end << "]" << std::endl;
	  out_ranges.push_back (newrange);
	  wire_ranges[i].first = wirerange.start;
	}
	// if we overlfow on the high-end
	if ( wirerange.end > wire_range.second ){
	  found = true;
	  btutil::WireRange_t newrange( wirerange.ch, wire_range.second, wirerange.end);
	  //std::cout << "\t adding wire-range @ chan "<< newrange.ch << " in range [" << newrange.start << ", " << newrange.end << "]" << std::endl;
	  out_ranges.push_back (newrange);
	  wire_ranges[i].second = wirerange.end;
	}
      }// for all wire-ranges already-found on this channel
      if (!found){
	// return the same exact wire-rang
	//std::cout << "\t this wire-range has not been added..." << std::endl;
	//std::cout << "\t adding wire-range @ chan "<< wirerange.ch << " in range [" << wirerange.start << ", " << wirerange.end << "]" << std::endl;
	out_ranges.push_back( wirerange );
	std::pair<double,double> tickrange = std::make_pair(wirerange.start,wirerange.end);
	std::vector< std::pair<double,double> > tickranges = {tickrange};
	_wiremap[ wirerange.ch ] = tickranges;
      }
    }// if there is some overlap
    
    return out_ranges;
  }

}
#endif
