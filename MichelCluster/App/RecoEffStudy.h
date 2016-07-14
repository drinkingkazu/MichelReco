/**
 * \file RecoEffStudy.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class RecoEffStudy
 *
 * @author david
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_RECOEFFSTUDY_H
#define LARLITE_RECOEFFSTUDY_H

#include "Analysis/ana_base.h"

#include "TTree.h"

#include <map>

namespace larlite {
  /**
     \class RecoEffStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class RecoEffStudy : public ana_base{
  
  public:

    /// Default constructor
    RecoEffStudy();

    /// Default destructor
    virtual ~RecoEffStudy(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetDebug(bool on) { _debug = on; }

  protected:

    // map linking track ID -> index in ev_track
    std::map<int,int> _trackIDMap;
    // vector containing mcshower index for Michel electrons
    std::vector<int> _michel_idx_v;
    // vector cointaing Michel electron start point [wire,tick]
    std::vector< std::pair<double,double> > _michel_start_v;

    TTree *_tree;
    double _mc_X, _mc_Y, _mc_Z, _mc_T;
    int    _mc_wire, _rc_wire;
    double _mc_tick, _rc_tick;

    double _mc_muon_E, _mc_muon_px, _mc_muon_py, _mc_muon_pz;
    double _mc_michel_E, _mc_michel_px, _mc_michel_py, _mc_michel_pz;

    double _rc_michel_E;

    double _3Ddot, _2Ddot;
    
    bool _debug;
    
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
