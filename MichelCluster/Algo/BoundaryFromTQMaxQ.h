/**
 * \file BoundaryFromTQMaxQ.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BoundaryFromTQMaxQ
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_BOUNDARYFROMTQMAXQ_H
#define MICHELCLUSTER_BOUNDARYFROMTQMAXQ_H

#include "Fmwk/BaseMichelAlgo.h"
#include <algorithm>

namespace michel {
  /**
     \class BoundaryFromTQMaxQ
  */
  class BoundaryFromTQMaxQ : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    BoundaryFromTQMaxQ()
      :
      _maxDistance(20)
    {_name="BoundaryFromTQMaxQ";}
    
    /// Default destructor
    ~BoundaryFromTQMaxQ(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// setter function for max distance between min in dQds
    /// and max in dQ
    void SetMaxDistancesTruncatedQMaxQ(size_t n) { _maxDistance = n; }
    
  private:

    /// this max distance is the window we should look in when we scan
    /// the true charge spectrum... we get the truncated max Q index
    /// then we go seach in true Q spectrum in _maxDistance window
    /// and find the maximum in that window. This becomes the michel
    /// boundary
    
    size_t _maxDistance;
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

