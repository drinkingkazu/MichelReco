#ifndef MICHELCLUSTER_BOUNDARYFROMTQMAXQ_CXX
#define MICHELCLUSTER_BOUNDARYFROMTQMAXQ_CXX

#include "BoundaryFromTQMaxQ.h"
#include "Fmwk/MichelException.h"
#include "ClusterVectorCalculator.h"
#include <cmath>
#include <sstream>

namespace michel {
  
  void BoundaryFromTQMaxQ::EventReset()
  {}
  
  bool BoundaryFromTQMaxQ::ProcessCluster(MichelCluster& cluster,
					  const std::vector<HitPt>& hits)
  { 

    ClusterVectorCalculator _clusterCalc;
    
    int candidate_loc     = _clusterCalc.find_max(cluster._t_mean_v); // Truncated max... hopefully not degenerate

    auto right = cluster._ordered_pts.size() - 1 - candidate_loc;
    auto left  = candidate_loc;
    
    int  iMin = 0;
    int  iMax = 0;    
    
    if(right >= _maxDistance) iMax  = _maxDistance   + candidate_loc;
    if(left  >= _maxDistance) iMin  = candidate_loc - _maxDistance;

    if(right < _maxDistance)  iMax  = cluster._hits.size() - 1;
    if(left  < _maxDistance)  iMin  = 0;

    // holder for hit with largest charge -> this will identify the location
    // of the hit that defines the michel boundary
    auto k   = 0.0;
    auto idx = 0;
    
    for(int w = iMin; w <= iMax; ++w) {
      auto c = cluster._hits[cluster._ordered_pts[w]]._q;
      // if this hit has more charge than any other
      if(c > k) { k = c; idx = w; }
    }
    
    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "Boundary set @ dS = " << cluster._s_v[idx];
      Print(msg::kINFO,this->Name(),ss.str());
    }

    cluster._boundary    = cluster._ordered_pts[idx];

    return true;
  }

}
#endif
