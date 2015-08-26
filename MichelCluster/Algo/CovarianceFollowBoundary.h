/**
 * \file CovarianceFollowBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CovarianceFollowBoundary
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_COVARAIANCEFOLLOWBOUNDARY_H
#define MICHELCLUSTER_COVARAIANCEFOLLOWBOUNDARY_H

#include "Fmwk/BaseAlgBoundary.h"
#include <algorithm>

namespace michel {
  /**
     \class CovarianceFollowBoundary
  */
  class CovarianceFollowBoundary : public BaseAlgBoundary {
    
  public:
    
    /// Default constructor
    CovarianceFollowBoundary() : _maxDistance(20),
				 _covariance_dip_cutoff(0.9),
				 _maxCovarianceAtStart(0.9)
    {}
    
    /// Default destructor
    ~CovarianceFollowBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(MichelCluster& cluster);

    /// setter function for max distance between min in dQds
    /// and max in dQ
    void SetMaxDistanceTruncatedPeaks(int n) { _maxDistance = n; }

    /// setter function for maxCovarianceAtStart
    void SetMaxCovarianceAtStart(double m)   { _maxCovarianceAtStart = m; }
    
    /// setter function for strength SINGLE dip in covariance must be
    /// 0.9 is default
    void SetCovarianceDipCutoff(double d)    { _covariance_dip_cutoff = d; }
    
  private:

    /// number of hits that sets max allowed distance between
    /// min in trunc. dQdx and max in trunc. dQ
    /// if the min and max are more than this number of hits
    /// apart -> do not create a michel
    int _maxDistance;

    /// single dip in covariance strenght, must be below this value to
    /// be considered a dip
    double _covariance_dip_cutoff;

    /// value that chosen start point covariance must be less than to
    /// be considered, stolen from David C. MatchBoundaries
    double _maxCovarianceAtStart;
    
    
    /// return sign of val (-1/1 and sometimes 0)
    int sign(double val);
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

