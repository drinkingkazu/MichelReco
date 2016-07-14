#ifndef LARLITE_RECOEFFSTUDY_CXX
#define LARLITE_RECOEFFSTUDY_CXX

#include "RecoEffStudy.h"

#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/cluster.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"

namespace larlite {

  RecoEffStudy::RecoEffStudy()
    : _tree(nullptr)
  {
    _debug = false;
    _name  = "RecoEffStudy";
    _fout  = 0;
  }

  bool RecoEffStudy::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","tree");
    _tree->Branch("_mc_X",&_mc_X,"mc_X/D");
    _tree->Branch("_mc_Y",&_mc_Y,"mc_Y/D");
    _tree->Branch("_mc_Z",&_mc_Z,"mc_Z/D");
    _tree->Branch("_mc_T",&_mc_T,"mc_T/D");
    _tree->Branch("_mc_wire",&_mc_wire,"mc_wire/I");
    _tree->Branch("_rc_wire",&_rc_wire,"rc_wire/I");
    _tree->Branch("_mc_tick",&_mc_tick,"mc_tick/D");
    _tree->Branch("_rc_tick",&_rc_tick,"rc_tick/D");

    _tree->Branch("_mc_muon_E",&_mc_muon_E,"mc_muon_E/D");
    _tree->Branch("_mc_muon_px",&_mc_muon_px,"mc_muon_px/D");
    _tree->Branch("_mc_muon_py",&_mc_muon_py,"mc_muon_py/D");
    _tree->Branch("_mc_muon_pz",&_mc_muon_pz,"mc_muon_pz/D");
    _tree->Branch("_mc_michel_E",&_mc_michel_E,"mc_michel_E/D");
    _tree->Branch("_mc_michel_px",&_mc_michel_px,"mc_michel_px/D");
    _tree->Branch("_mc_michel_py",&_mc_michel_py,"mc_michel_py/D");
    _tree->Branch("_mc_michel_pz",&_mc_michel_pz,"mc_michel_pz/D");

    _tree->Branch("_rc_michel_E",&_rc_michel_E,"rc_michel_E/D");

    _tree->Branch("_3Ddot",&_3Ddot,"3Ddot/D");
    _tree->Branch("_2Ddot",&_2Ddot,"2Ddot/D"); // w.r.t. collection plane

    return true;
  }
  
  bool RecoEffStudy::analyze(storage_manager* storage) {
    
    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    double w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    // time->cm (accounting for different operating voltages)
    double driftVel = larutil::LArProperties::GetME()->DriftVelocity(0.273,87); // [cm/us]
    // tick width in time
    double tickWidth = 0.5; // [us]
    double t2cm = tickWidth*driftVel;

    _trackIDMap.clear();
    _michel_start_v.clear();

    // load MCShower / MCTracks
    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto ev_mctrack  = storage->get_data<event_mctrack>("mcreco");

    // load Michel clusters
    auto ev_cluster  = storage->get_data<event_cluster>("michel");

    if (_debug)
      std::cout << "found " << ev_cluster->size() << " michels" << std::endl;
    
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
      
      // if no energy deposited in TPC by shower -> exit
      if (mcsh.DetProfile().E() < 1)
	continue;

      // get the parent nuon
      auto const& muon = ev_mctrack->at( _trackIDMap[ mcsh.MotherTrackID() ] );

      // if muon's size is < 2 points -> ignore...
      if (muon.size() < 2)
	continue;

      auto const& mu_end = muon.at( muon.size() - 1 );

      // if distance between muon and mcshower is to large -> ignore as muon decay
      //if ( sqrt ( pow(mu_end.X() - e_strt.X(),2) + pow(mu_end.Y() - e_strt.Y(),2) + pow(mu_end.Z() - e_strt.Z(),2) )  > 1 )
      //continue;

      // finally, we have a Michel electron!

      _mc_X = e_strt.X();
      _mc_Y = e_strt.Y();
      _mc_Z = e_strt.Z();
      _mc_T = e_strt.T();

      // project the Michel's start point onto the collection plane
      double start_w = _mc_Z / w2cm;
      double start_t = _mc_X / t2cm;
      // account for offset in trigger [T0]
      // get time in us and get distance w/ drift-velocity
      start_t += (_mc_T / 1000.) * 0.11 / t2cm;

      if (_debug)
	std::cout << "Found Michel starting @ [X,Z,T] -> [" << _mc_X << ", " << _mc_Z
		  << ", "<<  (double)_mc_T / 1000. << "]" << std::endl
		  << "Corresponding to [W,T] -> [" << start_w <<  ", " << start_t << "]" << std::endl
		  << "With w2cm = " << w2cm << " and t2cm " << t2cm << std::endl << std::endl;
      
      // save it's index
      _michel_idx_v.push_back( i );
      // save start point info
      _michel_start_v.push_back( std::make_pair(start_w,start_t) );
      
    }// for all mcshowers
    
    // load through Michel clusters and find the closest truth Michel to it
    for (size_t j=0; j < ev_cluster->size(); j++){
      
      auto const& clus = ev_cluster->at(j);

      _rc_wire = (int)clus.StartWire();
      _rc_tick = clus.StartTick();

      if (_debug)
	std::cout << "Michel cluster w/ [W,T] -> [" << _rc_wire
		  << ", " << _rc_tick << "]" << std::endl << std::endl;

      // find closest michel
      double d_min = 10000.;
      
      // loop through MC michels and find the one that best matches
      for (auto const& mc_michel : _michel_start_v){

	if ( (mc_michel.first - _rc_wire) * (mc_michel.first - _rc_wire)  < d_min){
	  d_min = (mc_michel.first - _rc_wire) * (mc_michel.first - _rc_wire);
	  _mc_wire = (int)mc_michel.first;
	  _mc_tick = mc_michel.second;
	}

      }// for all MC michels

      _tree->Fill();
      
    }// for all Michel clusters

    
    return true;
  }

  bool RecoEffStudy::finalize() {

    _tree->Write();

    return true;
  }

}
#endif
