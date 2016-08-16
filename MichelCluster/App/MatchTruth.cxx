#ifndef MATCHTRUTH_CXX
#define MATCHTRUTH_CXX

#include "MatchTruth.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"

namespace larlite{

  MatchTruth::MatchTruth()
  {
    _debug = true;
  }

  void MatchTruth::Match(const larlite::event_mcshower* ev_mcshower,
			 const larlite::event_mctrack*  ev_mctrack,
			 const larlite::event_cluster*  ev_cluster)
  {

    Reset();

    // use instances of LArUtil and GeometryUtilities
    // for (w,t) -> (cm, cm) conversion
    // wire->cm
    _w2cm = larutil::GeometryHelper::GetME()->WireToCm();
    _t2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    double efield   = larutil::LArProperties::GetME()->Efield(); // kV/cm
    double temp     = larutil::LArProperties::GetME()->Temperature(); // Kelvin
    _driftVel       = larutil::LArProperties::GetME()->DriftVelocity(efield,temp); // [cm/us]

    // build ID -> position map for MCTracks
    for (size_t i=0; i < ev_mctrack->size(); i++)
      _trackIDMap[ ev_mctrack->at(i).TrackID() ] = i;
    
    for (size_t i=0; i < ev_mcshower->size(); i++){
      
      auto const& mcsh = ev_mcshower->at(i);
      
      // only keep decays
      if ( (mcsh.Process() != "Decay") and (mcsh.Process() != "muMinusCaptureAtRest") )
	continue;

      if (mcsh.PdgCode() == 22)
	continue;
      
      // with parent a muon
      if ( (mcsh.MotherPdgCode() != 13) and (mcsh.MotherPdgCode() != -13) )
	continue;
      
      auto const& e_strt = mcsh.Start();

      if (e_strt.E() < 1)
	continue;
      
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

      if (_debug)
	std::cout << "Found a true michel w/ energy " << e_strt.E()
		  << " and process " << mcsh.Process()
		  << " parent muon trackID : " << mcsh.MotherTrackID() 
		  << std::endl;
    
      // project the Michel's start point onto the collection plane
      double start_w = e_strt.Z() / _w2cm;
      double start_t = e_strt.X() / _t2cm;
      // account for offset in trigger [T0]
      // get time in us and get distance w/ drift-velocity
      start_t += ( e_strt.T() / 1000.) * _driftVel / _t2cm + 800;

      // if start tick out of truncated waveform bounds -> don't include
      if ( (start_t < 0) or (start_t > 6300) )
	   continue;

      // save it's index
      if (_debug)
	std::cout << "found a michel @ ev_mcshower index " << i << std::endl;
      _mc_michel_idx_v.push_back( i );
      // save start point info
      _mc_michel_start_v.push_back( std::make_pair(start_w,start_t) );

      _muon_michel_idx_map[ _mc_michel_start_v.size() - 1 ] = std::make_pair( _trackIDMap[ mcsh.MotherTrackID() ], i );
      
    }// for all mcshowers

    std::cout << "what michels did we find?" << std::endl;
    for (size_t i=0; i < _mc_michel_idx_v.size(); i++)
      std::cout << "\tindex " << _mc_michel_idx_v[i] << std::endl;


    //  save information on reconstructed michels
    if (ev_cluster){
      for (size_t i=0; i < ev_cluster->size(); i++){
	auto const& clus = ev_cluster->at(i);
	_rc_michel_idx_v.push_back( i );
	_rc_michel_start_v.push_back( std::make_pair( (double)clus.StartWire(), clus.StartTick() ) );
	_rc_michel_elecQ_v.push_back( clus.SummedADC() );
      }
    }

    // we have selected all MC and RECO michel electrons -> now match
    // MatchMCtoReco();

    MatchRecotoMC();
    
    return;
  }

  // each MC michel is assigned a reco michel
  void MatchTruth::MatchMCtoReco()
  {

    _MCtoReco_match_v.clear();
    
    int    idx   = -1;
    double d_min = 10000.;

    for (size_t i=0; i < _mc_michel_start_v.size(); i++){

      auto const& mc_michel_start = _mc_michel_start_v[ i ];

      for (size_t j=0; j < _rc_michel_start_v.size(); j++){

	auto const& rc_michel_start = _rc_michel_start_v[ j ];

	double d = mc_michel_start.first - rc_michel_start.first;
	d *= d;

	if (d < d_min){
	  idx = j;
	  d_min = d;
	}// if closest

      }// for all reco michels

      // save this index
      _MCtoReco_match_v[ _mc_michel_idx_v [ i ] = _rc_michel_idx_v [ idx ] ];

    }// for all MC michels

    return;
  }


    // each MC michel is assigned a reco michel
  void MatchTruth::MatchRecotoMC()
  {

    _RecotoMC_match_v.clear();
    
    int    idx   = -1;
    double d_min = 10000.;

    for (size_t i=0; i < _rc_michel_start_v.size(); i++){

      auto const& rc_michel_start = _rc_michel_start_v[ i ];

      for (size_t j=0; j < _mc_michel_start_v.size(); j++){

	auto const& mc_michel_start = _mc_michel_start_v[ j ];

	double d = mc_michel_start.first - rc_michel_start.first;
	d *= d;

	if (_debug)
	  std::cout << " idx [i, j] -> [" << i << ", " << j << "] have distance " << d << std::endl;
	
	if (d < d_min){
	  if (_debug) std::cout << "\t closest..." << std::endl;
	  idx = j;
	  d_min = d;
	}// if closest

      }// for all reco michels

      if (idx == -1){
	if (_debug) std::cout << "no close matching michel found -> don't add" << std::endl;
	continue;
      }

      std::cout << "what michels did we find?" << std::endl;
      for (size_t k=0; k < _mc_michel_idx_v.size(); k++)
	std::cout << "\tindex " << _mc_michel_idx_v[k] << std::endl;
      
      // save this index
      if (_debug)
	std::cout << "matched pair : [i,j] -> [ " << i << ", " << idx 
		  << "] corresponds to indices [ " << _rc_michel_idx_v[i]
		  << ", " << _mc_michel_idx_v[idx] << "]" << std::endl;
      _RecotoMC_match_v[ _rc_michel_idx_v [ i ] ] = _mc_michel_idx_v [ idx ];

    }// for all MC michels
    
    return;
  }

  void MatchTruth::Reset()
  {

    _trackIDMap.clear();
    _mc_michel_idx_v.clear();
    _rc_michel_idx_v.clear();
    _mc_michel_start_v.clear();
    _rc_michel_start_v.clear();
    _rc_michel_elecQ_v.clear();
    _muon_michel_idx_map.clear();

    return;
  }

}

#endif
