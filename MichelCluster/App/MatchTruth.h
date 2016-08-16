/**
 * \file MatchTruth.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class MatchTruth
 *
 * @author david
 */

/** \addtogroup App

    @{*/
#ifndef MATCHTRUTH_H
#define MATCHTRUTH_H

#include <iostream>
#include <map>
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/cluster.h"

/**
   \class MatchTruth
   User defined class MatchTruth ... these comments are used to generate
   doxygen documentation!
 */

namespace larlite{
  
  class MatchTruth{

  public:
    
    /// Default constructor
    MatchTruth();
    
    /// Default destructor
    ~MatchTruth(){}

    void Match(const larlite::event_mcshower* ev_mcshower,
	       const larlite::event_mctrack*  ev_mctrack,
	       const larlite::event_cluster*  ev_cluster);

    void SetDebug(bool on) { _debug = on; }

    // grab matched pairs
    const std::map<int, int> GetRecotoMCMatch() { return _RecotoMC_match_v; }
    const std::map<int, int> GetMCtoRecoMatch() { return _MCtoReco_match_v; }

  private:

    void Reset();

    void MatchMCtoReco();
	
    void MatchRecotoMC();


    // map linking track ID -> index in ev_track
    std::map<int,int> _trackIDMap;
    // vector containing mcshower index for Michel electrons
    std::vector<int> _mc_michel_idx_v;
    // vector containing reconstructed michel index for Michel electrons
    std::vector<int> _rc_michel_idx_v;
    // vector cointaing Michel electron start point [wire,tick] from MC
    std::vector< std::pair<double,double> > _mc_michel_start_v;
    // from reco
    std::vector< std::pair<double,double> > _rc_michel_start_v;
    // vector containing ADC charge from electron-part only
    std::vector< double > _rc_michel_elecQ_v;
    // map linking _michel_start_v idx to [mu idx, michel idx];
    std::map<int, std::pair<int,int> > _muon_michel_idx_map;

    // match MC to reco
    // links index in event_mcshower to event_cluster
    std::map<int, int> _MCtoReco_match_v;
    // match reco to MC
    std::map<int, int> _RecotoMC_match_v;

    double _t2cm, _w2cm, _driftVel;

    bool _debug;
    
  };

}

#endif
/** @} */ // end of doxygen group 

