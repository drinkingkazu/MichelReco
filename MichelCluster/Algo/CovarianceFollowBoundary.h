/**
 * \file CovarianceFollowBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CovarianceFollowBoundary
 *
 * @author kazuhiro
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
    CovarianceFollowBoundary(){_maxDistance=20;}
    
    /// Default destructor
    ~CovarianceFollowBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(MichelCluster& cluster);

    /// setter function for max distance between min in dQds
    /// and max in dQ
    void SetMaxDistanceTruncatedPeaks(int n) { _maxDistance = n; }

    
  private:

    /// number of hits that sets max allowed distance between
    /// min in trunc. dQdx and max in trunc. dQ
    /// if the min and max are more than this number of hits
    /// apart -> do not create a michel
    int _maxDistance;

    
    int sign(double val);
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

