#ifndef LARLITE_RECOEFFSTUDY_CXX
#define LARLITE_RECOEFFSTUDY_CXX

#include "RecoEffStudy.h"


#include "DataFormat/cluster.h"
#include "DataFormat/trigger.h"
#include "DataFormat/hit.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"

namespace larlite {

  RecoEffStudy::RecoEffStudy()
    : _tree(nullptr)
  {
    _debug  = false;
    _name   = "RecoEffStudy";
    _all_mc = true;
    _fout   = 0;
    _distance = 3;
  }

  bool RecoEffStudy::initialize() {

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    _w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    _t2cm = larutil::GeometryHelper::GetME()->TimeToCm();



    if (_tree) delete _tree;

    _tree = new TTree("_tree","tree");
    _tree->Branch("_mc_X",&_mc_X,"mc_X/D");
    _tree->Branch("_mc_Y",&_mc_Y,"mc_Y/D");
    _tree->Branch("_mc_Z",&_mc_Z,"mc_Z/D");
    _tree->Branch("_mc_T",&_mc_T,"mc_T/D");
    _tree->Branch("_mc_wire",&_mc_wire,"mc_wire/D");
    _tree->Branch("_rc_wire",&_rc_wire,"rc_wire/D");
    _tree->Branch("_mc_tick",&_mc_tick,"mc_tick/D");
    _tree->Branch("_mc_tick_muon",&_mc_tick_muon,"mc_tick_muon/D");
    _tree->Branch("_rc_tick",&_rc_tick,"rc_tick/D");
    _tree->Branch("_rc_ADCq",&_rc_ADCq,"rc_ADCq/D");

    _tree->Branch("_mc_muon_E",&_mc_muon_E,"mc_muon_E/D");
    _tree->Branch("_mc_muon_px",&_mc_muon_px,"mc_muon_px/D");
    _tree->Branch("_mc_muon_py",&_mc_muon_py,"mc_muon_py/D");
    _tree->Branch("_mc_muon_pz",&_mc_muon_pz,"mc_muon_pz/D");

    _tree->Branch("_mc_michel_E",&_mc_michel_E,"mc_michel_E/D");
    _tree->Branch("_mc_michel_px",&_mc_michel_px,"mc_michel_px/D");
    _tree->Branch("_mc_michel_py",&_mc_michel_py,"mc_michel_py/D");
    _tree->Branch("_mc_michel_pz",&_mc_michel_pz,"mc_michel_pz/D");

    _tree->Branch("_mc_muon_decay_T",&_mc_muon_decay_T,"mc_muon_decay_T/D");
    _tree->Branch("_mc_michel_creation_T",&_mc_michel_creation_T,"mc_michel_creation_T/D");


    _tree->Branch("_rc_michel_E",&_rc_michel_E,"rc_michel_E/D");

    _tree->Branch("_trig_time",&_trig_time,"trig_time/D");

    _tree->Branch("_matched",&_matched,"matched/I");

    _tree->Branch("_3Ddot",&_3Ddot,"3Ddot/D");
    _tree->Branch("_2Ddot",&_2Ddot,"2Ddot/D"); // w.r.t. collection plane

    return true;
  }
  
  bool RecoEffStudy::analyze(storage_manager* storage) {

    double efield   = larutil::LArProperties::GetME()->Efield(); // kV/cm
    double temp     = larutil::LArProperties::GetME()->Temperature(); // Kelvin
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(efield,temp); // [cm/us]
    
    _trackIDMap.clear();
    _mc_michel_start_v.clear();
    _rc_michel_start_v.clear();
    _muon_michel_idx_map.clear();

    // load MCShower / MCTracks
    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto ev_mctrack  = storage->get_data<event_mctrack>("mcreco");

    // load Michel clusters
    auto ev_cluster  = storage->get_data<event_cluster>("michel");

    // load hits associated to Michels
    larlite::event_hit *ev_hit = nullptr;
    auto const& hit_ass_set = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());

    if (_debug)
      std::cout << std::endl << "found " << ev_cluster->size() << " michels" << std::endl;
    
    // build ID -> position map for MCTracks
    for (size_t i=0; i < ev_mctrack->size(); i++)
      _trackIDMap[ ev_mctrack->at(i).TrackID() ] = i;
    
    for (size_t i=0; i < ev_mcshower->size(); i++){
      
      auto const& mcsh = ev_mcshower->at(i);
      
      // only keep decays
      if ( (mcsh.Process() != "Decay") and (mcsh.Process() != "muMinusCaptureAtRest") )
	continue;

      // with parent a muon
      if ( (mcsh.MotherPdgCode() != 13) and (mcsh.MotherPdgCode() != -13) )
	continue;
      
      auto const& e_strt = mcsh.Start();
      
      // make sure shower starts in TPC
      if ( (e_strt.X() < 0) or (e_strt.X() > 256) or 
	   (e_strt.Y() < -116) or (e_strt.Y() > 116) or
	   (e_strt.Z() < 0) or (e_strt.Z() > 1036) )
	continue;
      
      // get the parent nuon
      auto const& muon = ev_mctrack->at( _trackIDMap[ mcsh.MotherTrackID() ] );

      // if muon's size is < 2 points -> ignore...
      if (muon.size() < 2)
	continue;

      // if muon does not decay at rest -> exit
      if ( muon.at(muon.size() - 2).E() > 106)
	continue;
    
      // project the Michel's start point onto the collection plane
      double start_w = e_strt.Z() / _w2cm;
      double start_t = e_strt.X() / _t2cm;
      // account for offset in trigger [T0]
      // get time in us and get distance w/ drift-velocity
      start_t += ( e_strt.T() / 1000.) * driftVel / _t2cm + 800;

      // if start tick out of truncated waveform bounds -> don't include
      if ( (start_t < 0) or (start_t > 6300) )
	   continue;

      if (_debug)
	std::cout << "Found Michel starting @ [X,Z,T] -> [" << _mc_X << ", " << _mc_Z
		  << ", "<<  (double)_mc_T / 1000. << "]" << std::endl
		  << "Corresponding to [W,T] -> [" << start_w <<  ", " << start_t << "]" << std::endl
		  << "With w2cm = " << _w2cm << " and t2cm " << _t2cm << std::endl << std::endl;
      
      // save it's index
      _michel_idx_v.push_back( i );
      // save start point info
      _mc_michel_start_v.push_back( std::make_pair(start_w,start_t) );

      _muon_michel_idx_map[ _mc_michel_start_v.size() - 1 ] = std::make_pair( _trackIDMap[ mcsh.MotherTrackID() ], i );
      
    }// for all mcshowers

    //  save information on reconstructed michels
    if (ev_cluster){
      for (auto const& clus : *ev_cluster){
	_rc_michel_start_v.push_back( std::make_pair( (double)clus.StartWire(), clus.StartTick() ) );
      }
    }

    // **************************
    // SAVE ENTRY PER TRUE MICHEL
    // **************************
    if (_all_mc){

      for (size_t i=0; i < _mc_michel_start_v.size(); i++){
	
	ResetTTree();
	
	// grab Michel MCShower
	auto const& michel_MCShower = ev_mcshower->at( _muon_michel_idx_map[ i ].second );
	auto const& muon_MCTrack    = ev_mctrack ->at( _muon_michel_idx_map[ i ].first  );
	
	if (muon_MCTrack.size() < 2)
	  continue;
	
	FillMuonInfo(muon_MCTrack);
	
	FillMichelInfo(michel_MCShower);
	
	FillDotProduct();
	
	// find best match from Reco'd michels
	auto const& matched = findBestMatch( _mc_michel_start_v[i], _rc_michel_start_v);
	
	_mc_tick = _mc_michel_start_v[i].second;
	_mc_wire = _mc_michel_start_v[i].first;
	
	if (matched.first == false){
	  _tree->Fill();
	  continue;
	}
	
	// if we found a good match, start filling info for RECO stuff
	auto const& matched_michel_cluster = ev_cluster->at( matched.second );
	_rc_wire = (double)matched_michel_cluster.StartWire();
	_rc_tick = matched_michel_cluster.StartTick();
	// grab hits associated to this michel
	auto const& hit_idx_v = hit_ass_set[ matched.second ];
	_rc_ADCq = 0.;
	for (auto const& hit_idx : hit_idx_v)
	  _rc_ADCq += ev_hit->at(hit_idx).Integral();
	_matched = 1;
	_tree->Fill();
      }// for all MC michels
      
    }// if we fill an entry per MC michel
    
    // **************************
    // SAVE ENTRY PER RECO MICHEL
    // **************************
    else{
      
      for (size_t i=0; i < _rc_michel_start_v.size(); i++){
	
	ResetTTree();
	
	// find best match from Reco'd michels
	auto const& matched = findBestMatch( _rc_michel_start_v[i], _mc_michel_start_v);

	_rc_tick = _rc_michel_start_v[i].second;
	_rc_wire = _rc_michel_start_v[i].first;

	if (matched.first == false){
	  _tree->Fill();
	  continue;
	}
	
	// grab Michel MCShower
	auto const& michel_MCShower = ev_mcshower->at( _muon_michel_idx_map[ matched.second ].second );
	auto const& muon_MCTrack    = ev_mctrack ->at( _muon_michel_idx_map[ matched.second ].first  );
	
	if (muon_MCTrack.size() < 2)
	  continue;
	
	FillMuonInfo(muon_MCTrack);
	
	FillMichelInfo(michel_MCShower);
	
	FillDotProduct();

	_mc_tick = _rc_michel_start_v[ matched.second ].second;
	_mc_wire = _rc_michel_start_v[ matched.second ].first;
	
	_matched = 1;
	_tree->Fill();
      }// for all MC michels
      
    }// if we fill an entry per MC michel
	
    return true;
  }

  bool RecoEffStudy::finalize() {

    _tree->Write();

    return true;
  }

  std::pair<bool, int> RecoEffStudy::findBestMatch(const std::pair<double,double>& coordinates,
						   const std::vector< std::pair<double,double> >& coord_v)
  {

    bool   found = false;
    int    idx   = -1;
    double d_min = 10000.;
    
    for (size_t i=0; i < coord_v.size(); i++){
      
      double d = coordinates.first - coord_v[i].first;
      d *= d;
      
      if (d < d_min){
	idx = i;
	d_min = d;
      }

    }// for all pairs of coordinates
    
    if (d_min < _distance*_distance)
      found = true;

    return std::make_pair(found, idx);
  }

  void RecoEffStudy::FillMuonInfo(const larlite::mctrack& muon)
  {

    _mc_muon_E  = muon.Start().E();
    _mc_muon_px = muon.Start().Px();
    _mc_muon_py = muon.Start().Py();
    _mc_muon_pz = muon.Start().Pz();


    
    double mag = sqrt( (_mc_muon_px * _mc_muon_px) +
		       (_mc_muon_py * _mc_muon_py) +
		       (_mc_muon_pz * _mc_muon_pz) );
		       
    _mc_muon_px /= mag;
    _mc_muon_py /= mag;
    _mc_muon_pz /= mag;

    auto const& mu_end = muon.at( muon.size() - 2 );
    
    _mc_muon_decay_T = mu_end.T();
    
    _mc_tick_muon = mu_end.X() / _t2cm + (mu_end.T() / 1000.) * ( 0.11 / _t2cm) ;

    return;
  }

  void RecoEffStudy::FillMichelInfo(const larlite::mcshower& michel)
  {

    auto const& e_strt = michel.Start();
    
    _mc_X = e_strt.X();
    _mc_Y = e_strt.Y();
    _mc_Z = e_strt.Z();
    _mc_T = e_strt.T();

    _mc_michel_E  = michel.Start().E();
    _mc_michel_px = michel.Start().Px();
    _mc_michel_py = michel.Start().Py();
    _mc_michel_pz = michel.Start().Pz();

    double mag = sqrt( (_mc_michel_px * _mc_michel_px) +
		       (_mc_michel_py * _mc_michel_py) +
		       (_mc_michel_pz * _mc_michel_pz) );
		       
    _mc_michel_px /= mag;
    _mc_michel_py /= mag;
    _mc_michel_pz /= mag;
    
    _mc_michel_creation_T = michel.DetProfile().T();

    return;
  }

  void RecoEffStudy::FillDotProduct()
  {

    _3Ddot = _mc_muon_px * _mc_michel_px + _mc_muon_py * _mc_michel_py  + _mc_muon_pz * _mc_michel_pz;
    
    // projected on XZ plane which is what is seen on collection plane
    _2Ddot = _mc_muon_px * _mc_michel_px + _mc_muon_pz * _mc_michel_pz;
    _2Ddot /= sqrt( _mc_muon_px * _mc_muon_px + _mc_muon_pz * _mc_muon_pz );
    _2Ddot /= sqrt( _mc_michel_px * _mc_michel_px + _mc_michel_pz * _mc_michel_pz );

  }

  void RecoEffStudy::ResetTTree()
  {

    _mc_X = _mc_Y = _mc_Z = _mc_T = kINVALID_DOUBLE;
    _mc_wire = _rc_wire = _mc_tick = _rc_tick = kINVALID_DOUBLE;
    _mc_tick_muon = kINVALID_DOUBLE;

    _rc_ADCq = 0;
    
    _trig_time = kINVALID_DOUBLE;
    
    _mc_muon_E = _mc_muon_px = _mc_muon_py = _mc_muon_pz = kINVALID_DOUBLE;
    _mc_michel_E = _mc_michel_px = _mc_michel_py = _mc_michel_pz = kINVALID_DOUBLE;
    
    _rc_michel_E = kINVALID_DOUBLE;
    
    _mc_muon_decay_T = _mc_michel_creation_T = kINVALID_DOUBLE;
    
    _3Ddot = _2Ddot = kINVALID_DOUBLE;

    _matched = 0;

  }
  
}
#endif
