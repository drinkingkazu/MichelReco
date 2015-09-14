/**
 * \file RequireBoundaryInLowCov.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RequireBoundaryInLowCov
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_COVARAIANCEFOLLOWBOUNDARY_H
#define MICHELCLUSTER_COVARAIANCEFOLLOWBOUNDARY_H

#include "Fmwk/BaseMichelAlgo.h"
#include <algorithm>

namespace michel {
  /**
     \class RequireBoundaryInLowCov
  */
  class RequireBoundaryInLowCov : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RequireBoundaryInLowCov()
      :
      _maxCovarianceAtStart(0.9)
    {_name="RequireBoundaryInLowCov";}
    
    /// Default destructor
    ~RequireBoundaryInLowCov(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);
    /// setter function for maxCovarianceAtStart
    void SetMaxCovarianceAtStart(double m)   { _maxCovarianceAtStart = m; }
    
  private:

    /// value that chosen start point covariance must be less than to
    /// be considered, stolen from David C. MatchBoundaries
    double _maxCovarianceAtStart;
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

