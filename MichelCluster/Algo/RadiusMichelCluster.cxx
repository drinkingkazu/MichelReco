#ifndef RADIUSMICHELCLUSTER_CXX
#define RADIUSMICHELCLUSTER_CXX

#include "RadiusMichelCluster.h"
#include <sstream>
namespace michel {

  void RadiusMichelCluster::EventReset()
  {}
  
  bool RadiusMichelCluster::ProcessCluster(MichelCluster& cluster,
					   const std::vector<HitPt>& hits) {

    if (hits.size() == 0) return false;
    
    auto& michel = cluster._michel;
    
    /// This michel was bogus when it came in, don't cluster further
    if (michel.size() == 0) return false;
    
    double radius = michel._length;
    
    /// Set radius to min rad of michel isn't long enough
    if (radius < _min_radius) radius = _min_radius;
    
    if (_verbosity <=  msg::kINFO) {
      std::stringstream ss;
      ss << "\tRadius clustering" << std::endl
	 << "\twe use a radius of: " << radius << std::endl
	 << "\tto cluster the remaining " << std::endl
	 << "\thits.size() : " << hits.size() << std::endl
	 << "\tinput";
      Print(msg::kINFO,__FUNCTION__,ss.str());
    }
    
    
    auto michel_start = michel._start;
    
    for (const auto& thishit : hits)
      if (michel_start.SqDist(thishit) <= radius * radius)
	if (thishit._q < _max_hit_charge)
	  michel.push_back(thishit);
    
    michel._charge = 0;
    
    for (const auto& michel_hit : michel)
      michel._charge += michel_hit._q;
   
    return true;
  }
  

  
}
#endif
