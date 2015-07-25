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
#include "BaseAlgIdentifier.h"
#include "MichelAnaBase.h"
#include <TFile.h>
#include <TStopwatch.h>
namespace michel {
  /**
     \class MichelReco
     User defined class MichelReco ... these comments are used to generate
     doxygen documentation!
  */
  class MichelReco{
    
  public:
    
    /// Default constructor
    MichelReco();
    
    /// Default destructor
    ~MichelReco(){}

    /// cluster register function
    void Append(const std::vector<michel::HitPt>& hit_v);

    /// Reco algorithm setter
    void SetAlgo(const AlgoType_t type, BaseMichelAlgo *algo);

    /// Ana algorithm adder
    void AddAna(MichelAnaBase* ana);

    /// Cluster configuration paramters
    void SetClusterConfig(size_t min_nhits, double d_cutoff)
    { _min_nhits = min_nhits; _d_cutoff = d_cutoff; }

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level) { _verbosity = level; }

    /// Initializer (before event loop)
    void Initialize();
    
    /// Executor (per event)
    void Process();

    /// Resetter (per event)
    void EventReset();

    /// Finalizer (after event loop)
    void Finalize(TFile *fout);

  protected:
    // MichelCluster configuration parameters
    double _d_cutoff;  ///< MichelCluster's cut-off distance for neighboring cluster
    size_t _min_nhits; ///< MichelCluster's min # hits to claim a cluster
    /// Verbosity
    msg::MSGLevel_t _verbosity;
    /// Algorithms to be executed
    std::vector< michel::BaseMichelAlgo* > _alg_v;
    /// Analysis to be executed
    std::vector< michel::MichelAnaBase* > _ana_v;
    /// Input clusters
    MichelClusterArray _input_v;
    /// Output clusters
    MichelClusterArray _output_v;
    //
    // Algorithms
    //
    BaseAlgMerger*     _alg_merge;    ///< Merging algorithm 
    BaseAlgBoundary*   _alg_boundary; ///< Michel/Muon boundary finder algorithm
    BaseAlgIdentifier* _alg_michel;   ///< Michel identification algorithm
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

