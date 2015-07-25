/**
 * \file MichelCluster.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class MichelCluster
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/
#ifndef MICHELCLUSTER_H
#define MICHELCLUSTER_H

//C++
#include <iostream>
#include <vector>
#include "MichelTypes.h"
#include "MichelConstants.h"
#include "Michel.h"
#include "HitPt.h"

namespace michel {
  /**
     \class MichelCluster
     User defined class MichelCluster ... these comments are used to generate
     doxygen documentation!
  */
  
  class MichelCluster {
    
  public:
    /// Default constructor
    MichelCluster(size_t min_nhits = 0,
		  double d_cutoff  = kMAX_DOUBLE );

    /// Default destructor
    virtual ~MichelCluster(){}

    /// Hit list setter
    void SetHits(const std::vector<michel::HitPt>& hits);

    /// Verbosity level
    void SetVerbosity(MSGLevel_t level)
    { _verbosity = level; }

    //
    // Operator attributes
    //
    // Addition
    MichelCluster operator+(const MichelCluster& other);
    
    //
    // Function attributes
    //
    /// Construct ordered index vector for near-by neighbor hits
    void OrderPoints(size_t start_index,
		     std::vector<size_t>& ordered_index_v,
		     std::vector<double>& ds_v,
		     std::vector<double>& s_v);

    /// Check if two clusters are touching 
    bool Touching (const MichelCluster& rhs,
		   const double min_dist);

    /// Find the closest hit to the reference from the cluster hit list
    const HitPt& ClosestHit(const HitPt& ref);

    /// Dumps information about this cluster
    void Dump();

    //
    // Data attributes
    //
    // Basic hit-based parameters
    std::vector<HitPt> _hits; ///< List of hits
    HitPt _start;             ///< Start point
    HitPt _end;               ///< End point

    // Ordered per-hit or in-between-hit quantities
    std::vector<HitIdx_t> _ordered_pts; ///< Distance ordered points
    std::vector<double>   _ds_v;        ///< Distance in-between neighboring hits (follows _ordered_pts) 
    std::vector<double>   _s_v;         ///< Cumulative neighboring distance between hits

    //
    // Attributes to be filled by external algorithms
    //
    HitIdx_t _boundary;   ///< Michel/Muon boundary point (not necessarily michel start exactly)
    Michel _michel;       ///< Michel    
    std::vector<double>   _chi2_v;      ///< Local linear chi2 fit 
    std::vector<double>   _t_mean_v;    ///< Truncated mean
    std::vector<double>   _t_dqds_v;    ///< Truncated dqds

  protected:
    //
    // Configurations
    //
    MSGLevel_t _verbosity; ///< Verbosity level
    size_t _min_nhits;     ///< Minimum # of hits
    double _d_cutoff;      ///< Distance cut off value

  };

  /// A set of MichelCluster
  typedef std::vector<michel::MichelCluster> MichelClusterArray;
}
#endif
/** @} */ // end of doxygen group 

