/**
 * \file RequireCloseTruncatedPeaks.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RequireCloseTruncatedPeaks
 *
 * @author vic
 */

/** \addtogroup MichelCluster
    
    @{*/
#ifndef MICHELCLUSTER_REQUIRECLOSETRUNCATEDPEAKS_H
#define MICHELCLUSTER_REQUIRECLOSETRUNCATEDPEAKS_H

#include "Fmwk/BaseMichelAlgo.h"
#include <algorithm>

namespace michel {
  /**
     \class RequireCloseTruncatedPeaks
  */
  class RequireCloseTruncatedPeaks : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RequireCloseTruncatedPeaks()
      :
      _maxDistance(20)
    {_name="RequireCloseTruncatedPeaks";}
    
    /// Default destructor
    ~RequireCloseTruncatedPeaks(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// setter function for max distance between min in dQds
    /// and max in dQ
    void SetMaxDistanceTruncatedPeaks(int n) { _maxDistance = n; }
    
  private:

    /// number of hits that sets max allowed distance between
    /// min in trunc. dQdx and max in trunc. dQ
    /// if the min and max are more than this number of hits
    /// apart -> do not create a michel
    int _maxDistance;
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

