/**
 * \file MichelMCStudy.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class MichelMCStudy
 *
 * @author davidc1
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_MICHELMCSTUDY_H
#define LARLITE_MICHELMCSTUDY_H

#include "Analysis/ana_base.h"
//Backtracker
#include "MCComp/MCMatchAlg.h"
// matching
#include "MatchTruth.h"

namespace larlite {
  /**
     \class MichelMCStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class MichelMCStudy : public ana_base{
  
  public:

    /// Default constructor
    MichelMCStudy();

    /// Default destructor
    virtual ~MichelMCStudy(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// set debugger for mcq
    void SetDebugMCQ(bool on) { _debug_mcq = on; }

    void SetFillHitTree(bool on) { _fill_hit_tree = on; }

  protected:

    MatchTruth _MatchTruth;

    void Reset();

    /// boolean -> debug mc charge collection info
    bool _debug_mcq;

    /// boool -> fill hit by hit tree
    bool _fill_hit_tree;

    // btalg instance for the entire Michel shower
    ::btutil::MCMatchAlg _BTAlgShower;
    // btalg instance for only the michel electron particle
    ::btutil::MCMatchAlg _BTAlgPart;

    // map to inform of wire-coverage from hit ranges
    std::map< unsigned int, std::vector< std::pair<double,double> > > _wiremap;

    // check if wire range was already used
    std::vector<btutil::WireRange_t> getUnUsedWireRange(const btutil::WireRange_t& wirerange);

    // stuff for MC study
    // Tree to store the amplitude of all hits
    TTree* _MIP_tree;
    double _hit_charge;

    // mc tree
    TTree* _mc_tree;
    // MC information for Michel electron
    double _mc_energy;
    double _mc_x, _mc_y, _mc_z, _mc_px, _mc_py, _mc_pz;
    double _mc_yplane_angle;
    double _mu_energy;
    double _mu_x, _mu_y, _mu_z, _mu_px, _mu_py, _mu_pz;
    double _mu_yplane_angle;
    // muon-michel dot-product
    double _mu_michel_3dangle;
    // reconstruted energy information
    double _reco_energy;
    std::vector<double> _michel_hit_fracReco;
    std::vector<double> _michel_hit_QtotReco;
    std::vector<int>    _michel_hit_idxReco;
    std::vector<double> _michel_hit_fracMC;
    std::vector<double> _michel_hit_QtotMC;
    double _QMichelMC;            // michel charge from MCShower Charge() function
    double _QMichelRecoSimch_all; // sum of Q from hits that have been reco'd as michel (entire hit contribution)
    double _QMichelRecoSimch_shr; // sum of Q from hits that have been reco'd as michel (only component from Michel IDEs)
    double _QMichelShowerMCSimch_all; // sum of simch contribution from hits to full MCShower (entire hit contribution)
    double _QMichelShowerMCSimch_shr; // sum of simch contribution from hits to full MCShower (only component from Michel IDEs)
    double _QMichelPartMCSimch_all;   // sum of simch contribution from hits to michel part only (entire hit contribution)
    double _QMichelPartMCSimch_shr;   // sum of simch contribution from hits to michel part only (only component from Michel IDEs)
    double _f_RecoHitsQ_fromMichelSimch; // fraction of charge collected in reco'd michel hits actually from michel
    int _run;
    int _subrun;
    int _event;

    // mc hit tree
    TTree* _mc_hit_tree;
    double _hit_integral;
    double _hit_mc_q;
    
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
