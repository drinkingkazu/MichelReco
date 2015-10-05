/**
 * \file MuonClusterTagger.h
 *
 * \ingroup MichelOpTagger
 * 
 * \brief Class def header for a class MuonClusterTagger
 *
 * @author david
 */

/** \addtogroup MichelOpTagger

    @{*/

#ifndef LARLITE_MUONCLUSTERTAGGER_H
#define LARLITE_MUONCLUSTERTAGGER_H

#include "Analysis/ana_base.h"
#include "OpT0Finder/Base/FlashMatchManager.h"
#include "FindFlashMichel.h"
#include <TTree.h>

namespace larlite {
  /**
     \class MuonClusterTagger
     User custom analysis class made by SHELL_USER_NAME
   */
  class MuonClusterTagger : public ana_base{
  
  public:

    /// Default constructor
    MuonClusterTagger();

    /// Default destructor
    virtual ~MuonClusterTagger(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    ::flashana::FlashMatchManager& Manager() { return _mgr;}

    void UseMC(bool on) { _use_mc = on; }

    void UseY(bool on) { _use_y = on; }

    void SetClusterProducer(std::string s) { _clusProducer = s; }
    
    void SetEfield(double efield) { _Efield = efield; }

    /// add michel finding algorithm
    void AddMichelFinder(larlite::FindFlashMichel algo) { _FindFlashMichel = algo; }

    /// verbosity flag setter
    void SetVerbose(bool on) { _verbose = on; _FindFlashMichel.SetVerbose(on); }

  protected:

    /// verbosity flag
    bool _verbose;

    std::string _clusProducer;
    
    double _Efield; // kV/cm

    double _w2cm, _t2cm;

    ::flashana::FlashMatchManager _mgr;

    bool _use_mc;

    bool _use_y;
    
    TTree* _tree;
    double _npe;
    double _flash_x;
    double _flash_y;
    double _flash_z;
    double _tpc_x;
    double _tpc_y;
    double _tpc_z;
    double _score;
    double _flash_time;
    // _t_drift : time the track actually drifted in the TPC
    // this is = ( RO time ) - ( matched PMT time )
    double _t_drift;
    double _mc_time;
    double _mc_x;
    double _mc_y;
    double _mc_z;
    double _mc_dx;

    // tree containing muon->michel optical match info
    TTree* _match_tree;
    double _tpc_flash_dz; // distance in z between matched tpc track and pmt flash
    double _dy, _dz;
    double _muon_t, _michel_t;
    double _muon_pe, _michel_pe;

    // FindOpMichel algorithm
    FindFlashMichel _FindFlashMichel;
    
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
