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
    double _mc_energy;
    double _reco_energy;
    std::vector<double> _michel_hit_frac;
    std::vector<double> _michel_hit_Qtot;
    double _QMichelMC; // michel charge from MCShower Charge() function
    double _QMichel; // sum of Q from hits > 10% michel in EDep (only Michel contribution to charge)
    double _QMichelTot; // sum of Q from all hits in EDep (only Michel contribution to charge)
    double _QMichelReco; // sum of Q from hits that have been reco'd as michel
    double _totQHits; // sum of Q from hits > 10% michel in EDep
    double _lifetimeCorr; // factor by which charge needs to be corrected to get an accurate lifetime
    double _mc_x; // x-position where michel starts

    /// boolean to select if to save Michel Clusters or not
    bool _save_clusters;

    /// electric field strength [ kV/cm]
    double _Efield;

    /// minimum cluster size for an input cluster to be considered
    int _minClusSize;

  private:
    ::btutil::MCMatchAlg _BTAlg;

    // text file where to dump input clsuter information
    std::ofstream _out_txt_file;

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
