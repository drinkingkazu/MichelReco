/**
 * \file MichelRecoDriver.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class MichelRecoDriver
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/

#ifndef LARLITE_MICHELRECODRIVER_H
#define LARLITE_MICHELRECODRIVER_H

#include "Analysis/ana_base.h"
#include "Fmwk/MichelRecoManager.h"
#include "Fmwk/MichelTypes.h"
#include <TStopwatch.h>
#include <TTree.h>
#include <fstream>
//Backtracker
#include "MCComp/MCMatchAlg.h"

namespace larlite {

  /**
     \class MichelRecoDriver
   */
  class MichelRecoDriver : public ana_base{
  
  public:

    /// Default constructor
    MichelRecoDriver();

    /// Default destructor
    virtual ~MichelRecoDriver(){}

    /** IMPLEMENT in MichelRecoDriver.cc!
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MichelRecoDriver.cc! 
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MichelRecoDriver.cc! 
    */
    virtual bool finalize();

    /// Reco manager getter
    ::michel::MichelRecoManager& GetManager() { return _mgr; }

    /**
     * @brief Set boolean to decide if to save michels in output cluster
     */
    void saveMichelClusters(bool on) { _save_clusters = on; }
    
    /// Producer setter
    void SetClusterProducer(const std::string& name) { _producer = name; }

    /// Set electric field strength (kV/cm)
    void SetEField(double E) { _Efield = E; }

    /// set if to use MC info or not
    void SetUseMC(bool on) { _use_mc = on; }

    /// Specifier for a specific plane reconstruction
    void SetPlane(unsigned int p) { _reco_plane.resize(p+1,false); _reco_plane[p]=true; }

    /// Set the minimum number of hits for a cluster to
    /// be considered as a candidate MichelCluster
    void SetMinClusSize(int n) { _minClusSize = n; }

    /// set debugger for mcq
    void SetDebugMCQ(bool on) { _debug_mcq = on; }

  protected:

    /// Stopwatch for event time profiling
    TStopwatch _event_watch; ///< For profiling
    double _event_time;
    size_t _event_ctr;
    
    /// Reco manager
    ::michel::MichelRecoManager _mgr;

    /// Input cluster producer name string
    std::string _producer;

    /// boolean. Use MC info?
    bool _use_mc;
    /// boolean -> debug mc charge collection info
    bool _debug_mcq;
    
    /// Option to set specific-plane-only reco
    std::vector<bool> _reco_plane;

    /// Output analysis TTree ptr
    TTree* _hit_tree;
    /// info for tree
    std::vector<double> _q_v;
    std::vector<double> _w_v;
    std::vector<double> _t_v;
    std::vector<double> _p_v;
    
    int _run;
    int _subrun;
    int _event;

    std::vector<michel::EventID> _events_info;

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

    // mc hit tree
    TTree* _mc_hit_tree;
    double _hit_integral;
    double _hit_mc_q;

    /// boolean to select if to save Michel Clusters or not
    bool _save_clusters;

    /// electric field strength [kV/cm]
    double _Efield;

    /// minimum cluster size for an input cluster to be considered
    int _minClusSize;

  private:
    // btalg instance for the entire Michel shower
    ::btutil::MCMatchAlg _BTAlgShower;
    // btalg instance for only the michel electron particle
    ::btutil::MCMatchAlg _BTAlgPart;

    // text file where to dump input clsuter information
    std::ofstream _out_txt_file;

    // map to inform of wire-coverage from hit ranges
    std::map< unsigned int, std::vector< std::pair<double,double> > > _wiremap;

    // check if wire range was already used
    std::vector<btutil::WireRange_t> getUnUsedWireRange(const btutil::WireRange_t& wirerange);

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
