/**
 * \file MichelRecoManager.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class MichelRecoManager
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_MICHELRECOMANAGER_H
#define MICHELCLUSTER_MICHELRECOMANAGER_H

#include <iostream>
#include "MichelCluster.h"
#include "BaseAlgMerger.h"
#include "BaseMichelAlgo.h"
#include "MichelAnaBase.h"
#include "ColorPrint.h"
#include <TFile.h>
#include <TStopwatch.h>
#include <math.h>


namespace michel {
/**
   \class MichelRecoManager
   A class that handles event-wise michel electron identification and reconstruction. \n
   It takes a list of clusters as an input where each cluster is represented as       \n
   a collection of michel::HitPt. Internally each cluster is converted into another   \n
   cluster representation michel::MichelCluster. Then a chain of following algorithms \n
   are applied to reconstruct michel electron(s):                                     \n
                                                                                      \n
   1) Input cluster merging                                                           \n
   2) Muon/Michel bounary finding                                                     \n
   3) Michel Mis ID Filtering                                                         \n
   3) Michel cluster ID                                                               \n
   4) Michel re-clustering                                                            \n
                                                                                      \n
   in a consecutive order.
*/
class MichelRecoManager : public ColorPrint {

public:

    /// Default constructor
    MichelRecoManager();

    /// Default destructor
    ~MichelRecoManager() {}

    //
    // Configuration functions
    //

    /// Add Merging ALgo
    void AddMergingAlgo(BaseAlgMerger *algo);

    /// Reco algorithm setter
    void AddAlgo(BaseMichelAlgo *algo);

    /// Ana algorithm adder
    void AddAna(MichelAnaBase* ana);

    /// Cluster configuration paramters
    void SetClusterConfig(size_t min_nhits, double d_cutoff)
    { _min_nhits = min_nhits; _d_cutoff = d_cutoff; }

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

    /// Setter for event info
    void SetEventInfo(EventID id) { _id = id; }

    /// Setter for distance cutoff when making internal track
    void SetDCutoff(double d)     { _d_cutoff = d; }

    /// Setter for setting minimum number of hits to create MichelCluster
    void SetMinNHits(int n)       { _min_nhits = n; }

    /// Setter for a debug mode
    void SetDebug(bool doit = true) { _debug = doit; }

    /// algorithm verbosity setter
    //void SetAlgoVerbosity(msg::MSGLevel_t verbose) { _algo_verbosity = verbose; }

    /// Getter for input MichelClusterArray
    MichelClusterArray GetMergedClusters()
    { return _merged_v; }

    /// Get algorithm vector
    std::vector< michel::BaseMichelAlgo* > GetAlgoVector()
      { return _alg_v; }    

    /// get hit vector
    const std::vector< ::michel::HitPt > GetHitVector()
    { return _all_hit_v; }

    /// Getter for output MichelClusterArray
    const MichelClusterArray& GetResult()
    { return _output_v; }

protected:

    // MichelCluster configuration parameters
    ///< MichelCluster's cut-off distance for neighboring cluster
    ///< MAXIMUM distance. Set this to a huge number for very generous ordering
    ///< Set this to a very small number for very conservative ordering.
    double _d_cutoff;

    // algorithm's verbosity
    //msg::MSGLevel_t _algo_verbosity;

    size_t _min_nhits; ///< MichelCluster's min # hits to claim a cluster
    /// Debug mode
    bool _debug;
    /// Input clusters
    MichelClusterArray _input_v;
    /// Merged clusters
    MichelClusterArray _merged_v;
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

    /// Event information
    EventID _id;

};
}

#endif
/** @} */ // end of doxygen group

