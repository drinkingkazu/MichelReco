/**
 * \file MatchBoundaries.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class MatchBoundaries
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_MATCHBOUNDARIES_H
#define MICHELCLUSTER_MATCHBOUNDARIES_H

#include "Fmwk/BaseAlgBoundary.h"
#include <algorithm>

namespace michel {
  /**
     \class MatchBoundaries
  */
  class MatchBoundaries : public BaseAlgBoundary {
    
  public:
    
    /// Default constructor
    MatchBoundaries(){_maxDistance=20;_maxCovarianceAtStart=0.9;}
    
    /// Default destructor
    ~MatchBoundaries(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(MichelCluster& cluster);

    /// setter function for max distance between min in dQds
    /// and max in dQ
    void SetMaxDistanceTruncatedPeaks(int n) { _maxDistance = n; }

    /// set the max covariance in the neightborhood of the start point
    void SetMaxCovarianceAtStart(double c) { _maxCovarianceAtStart = c; }

    
  private:

    /// number of hits that sets max allowed distance between
    /// min in trunc. dQdx and max in trunc. dQ
    /// if the min and max are more than this number of hits
    /// apart -> do not create a michel
    int _maxDistance;

    /// Maximum allowed covariance value for hits next to
    /// the start point
    double _maxCovarianceAtStart;
    
  };
}

#endif
/** @} */ // end of doxygen group 

