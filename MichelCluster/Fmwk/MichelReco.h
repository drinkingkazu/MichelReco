/**
 * \file MichelReco.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class MichelReco
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_MICHELRECO_H
#define MICHELCLUSTER_MICHELRECO_H

#include <iostream>
#include "MichelCluster.h"
#include "BaseAlgMerger.h"
#include "BaseAlgBoundary.h"
#include "BaseAlgMichelID.h"
#include "BaseAlgMichelCluster.h"
#include "MichelAnaBase.h"
#include <TFile.h>
#include <TStopwatch.h>
namespace michel {
  /**
     \class MichelReco
     A class that handles event-wise michel electron identification and reconstruction. \n
     It takes a list of clusters as an input where each cluster is represented as       \n
     a collection of michel::HitPt. Internally each cluster is converted into another   \n
     cluster representation michel::MichelCluster. Then a chain of following algorithms \n
     are applied to reconstruct michel electron(s):                                     \n
                                                                                        \n
     1) Input cluster merging                                                           \n
     2) Muon/Michel bounary finding                                                     \n
     3) Michel cluster ID                                                               \n
     4) Michel re-clustering                                                            \n
                                                                                        \n
     in a consecutive order.
  */
  class MichelReco{
    
  public:
    
    /// Default constructor
    MichelReco();
    
    /// Default destructor
    ~MichelReco(){}

    //
    // Configuration functions
    //
    /// Reco algorithm setter
    void SetAlgo(const AlgoType_t type, BaseMichelAlgo *algo);

    /// Ana algorithm adder
    void AddAna(MichelAnaBase* ana);

    /// Cluster configuration paramters
    void SetClusterConfig(size_t min_nhits, double d_cutoff)
    { _min_nhits = min_nhits; _d_cutoff = d_cutoff; }

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level) { _verbosity = level; }

    //
    // Data register functions
    //
    /// cluster register function
    void Append(const std::vector<michel::HitPt>& hit_v);
    /// all-hit register function
    void RegisterAllHits(const std::vector<michel::HitPt>& all_hit_v);
#ifndef __CINT__
    /// cluster register function w/ std::move
    void Append(std::vector<michel::HitPt>&& hit_v);
    /// all-hit register function w/ std::move
    void RegisterAllHits(std::vector<michel::HitPt>&& all_hit_v);
#endif
    
    //
    // Driver functions (change state)
    //
    /// Initializer (before event loop)
    void Initialize();
    
    /// Executor (per event)
    void Process();

    /// Resetter (per event)
    void EventReset();

    /// Finalizer (after event loop)
    void Finalize(TFile *fout);

    /// Getter for MichelClusterArray
    const MichelClusterArray& GetResult()
    { return _output_v; }

  protected:
    // MichelCluster configuration parameters
    double _d_cutoff;  ///< MichelCluster's cut-off distance for neighboring cluster
    size_t _min_nhits; ///< MichelCluster's min # hits to claim a cluster
    /// Verbosity
    msg::MSGLevel_t _verbosity;
    /// Input clusters
    MichelClusterArray _input_v;
    /// Output clusters
    MichelClusterArray _output_v;
    /// "ALL" hit list
    std::vector< ::michel::HitPt > _all_hit_v;
    /// Used hit marker for "ALL" hit list
    std::vector< bool > _used_hit_marker_v;
    //
    // Algorithms
    //
    BaseAlgMerger*        _alg_merge;          ///< Merging algorithm 
    BaseAlgBoundary*      _alg_boundary;       ///< Michel/Muon boundary finder algorithm
    BaseAlgMichelID*      _alg_michel_id;      ///< Michel identification algorithm
    BaseAlgMichelCluster* _alg_michel_cluster; ///< Michel re-clustering algorithm
    /// Algorithms to be executed
    std::vector< michel::BaseMichelAlgo* > _alg_v;
    /// Analysis to be executed
    std::vector< michel::MichelAnaBase* >  _ana_v;

    //
    // Time profilers
    //
    TStopwatch _watch; ///< For profiling
    std::vector<double> _alg_time_v; ///< Overall time for processing
    std::vector<size_t> _alg_ctr_v;  ///< Overall number of clusters processed by algo

  };
}

#endif
/** @} */ // end of doxygen group 

