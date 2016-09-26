/**
 * \file MichelTruthStudies.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class MichelTruthStudies
 *
 * @author dcaratelli
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_MICHELTRUTHSTUDIES_H
#define LARLITE_MICHELTRUTHSTUDIES_H

#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

#include "Analysis/ana_base.h"
#include "TTree.h"

#include <map>

namespace larlite {
  /**
     \class MichelTruthStudies
     User custom analysis class made by SHELL_USER_NAME
   */
  class MichelTruthStudies : public ana_base{
  
  public:

    /// Default constructor
    MichelTruthStudies();

    /// Default destructor
    virtual ~MichelTruthStudies(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    // map linking track ID -> index in ev_track
    std::map<int,int> _trackIDMap;

    TTree *_tree_mc;

    int _event, _subrun, _run;

    double _mc_X, _mc_Y, _mc_Z, _mc_T;

    double _trig_time;

    // mc information
    int _pdg, _parent_pdg;
    double _parent_end_E;
    std::string _process;
    
    int _mu_trkID;

    double _mc_muon_E, _mc_muon_px, _mc_muon_py, _mc_muon_pz;
    double _mc_michel_px, _mc_michel_py, _mc_michel_pz;
    double _mc_michel_E_edep, _mc_michel_E_true, _mc_michel_E_elec;

    double _rc_michel_E;

    double _mc_michel_creation_T;

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
