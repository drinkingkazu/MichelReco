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

#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

#include "MatchTruth.h"

#include "TTree.h"

#include <map>
#include <climits>
#include <limits>

namespace larlite {
  /**
     \class RecoEffStudy
     User custom analysis class made by SHELL_USER_NAME
   */

  static const double kINVALID_DOUBLE = std::numeric_limits<double>::max();

  static const double kMAX_DOUBLE     = std::numeric_limits<double>::max();
  
  static const double kMIN_DOUBLE     = std::numeric_limits<double>::min();

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
    
    void SetGoodMatchDistance(double d) { _distance = d; }
    
    void SetMCShowerProducer(std::string s) { _mcshower_producer = s; }

  protected:

    MatchTruth _MatchTruth;

    std::string _mcshower_producer;

    // function to get true best matched pair
    std::pair<bool, int> findBestMatch(const std::pair<double,double>& coordinates,
				       const std::vector< std::pair<double,double> >& coord_v);

    // fill muon information for TTree
    void FillMuonInfo(const larlite::mctrack& muon);

    // fill michel information for TTree
    void FillMichelInfo(const larlite::mcshower& michel);

    // fill dot products
    void FillDotProduct();

    double _w2cm, _t2cm;

    // map linking track ID -> index in ev_track
    std::map<int,int> _trackIDMap;
    // vector containing mcshower index for Michel electrons
    std::vector<int> _michel_idx_v;
    // vector cointaing Michel electron start point [wire,tick] from MC
    std::vector< std::pair<double,double> > _mc_michel_start_v;
    // from reco
    std::vector< std::pair<double,double> > _rc_michel_start_v;
    // vector containing ADC charge from electron-part only
    std::vector< double > _rc_michel_elecQ_v;
    // vector containing ADC charge for various photons
    std::vector< std::vector< double> > _rc_michel_photonQ_v;
    // vector containing X/Z coordinate of electron segments
    std::vector< std::vector< double > > _rc_elec_w_coord_v;
    std::vector< std::vector< double > > _rc_elec_t_coord_v;
    // map linking _michel_start_v idx to [mu idx, michel idx];
    std::map<int, std::pair<int,int> > _muon_michel_idx_map;

    TTree *_tree_mc;
    TTree *_tree_rc;

    int _event, _subrun, _run;

    double _mc_X, _mc_Y, _mc_Z, _mc_T;
    double _mc_wire, _rc_wire;
    double _mc_tick, _rc_tick;
    double _mc_tick_muon;
    double _rc_ADCq_elec, _rc_ADCq_tot;

    std::vector<double> _photon_q_v;
    std::vector<double> _electron_w_v;
    std::vector<double> _electron_t_v;

    double _trig_time;

    // mc information
    int _pdg, _parent_pdg;
    double _parent_end_E;
    std::string _process;
    

    double _mc_muon_E, _mc_muon_px, _mc_muon_py, _mc_muon_pz;
    double _mc_michel_px, _mc_michel_py, _mc_michel_pz;
    double _mc_michel_E_edep, _mc_michel_E_true, _mc_michel_E_elec;

    double _rc_michel_E;

    double _mc_muon_decay_T, _mc_michel_creation_T;

    double _3Ddot, _2Ddot;

    int    _matched;
    
    // reset all TTree variables
    void ResetTTree();

    // debug mode?
    bool _debug;
    
    // distance required for a good MC-RECO match
    double _distance;

    // drift velocity
    double _driftVel;
    
    
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
