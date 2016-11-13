#ifndef MICHELCLUSTER_MAXQBOUNDARY_CXX
#define MICHELCLUSTER_MAXQBOUNDARY_CXX

#include "MaxQBoundary.h"

namespace michel {

  void MaxQBoundary::EventReset()
  {}
    
  bool MaxQBoundary::ProcessCluster(MichelCluster& cluster,
				    const std::vector<HitPt>& hits)
  { 

    if (hits.size() == 0) return false;
    
    std::vector<double> charges;
    charges.reserve(cluster._ordered_pts.size());
    
    for(const auto& o: cluster._ordered_pts)
      charges.push_back(cluster._hits[o]._q);
    
    auto idx = find_max(charges);

    cluster._boundary = cluster._ordered_pts[idx]; 
    return true;
  }
  
  
  size_t MaxQBoundary::find_max(const std::vector<double>& data) {
    
    auto the_max = double{0.0};
    size_t idx = kINVALID_SIZE;
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] > the_max) {
	the_max = data[i]; idx = i;
      }
    }
    
    return idx;
  }

}
#endif
