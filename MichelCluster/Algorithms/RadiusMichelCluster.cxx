#ifndef RADIUSMICHELCLUSTER_CXX
#define RADIUSMICHELCLUSTER_CXX

#include "RadiusMichelCluster.h"

namespace michel {

  void RadiusMichelCluster::EventReset()
  {}
  
  void RadiusMichelCluster::Cluster(Michel& michel,
				    const std::vector<HitPt>& hits){
    

    
    //This michel was bogus when it came in, don't cluster further
    if (michel.size() == 0) return;

    
    double min_radius = 10; //hardcoded for now
    double radius = michel._length;
    
    //set to min rad
    if (radius < min_radius) radius = min_radius;

    if(_verbosity <=  msg::kINFO) {
      std::cout << "\n\n\tRadius clustering"
		<< "\n\twe use a radius of: " << radius
		<< "\n\tto cluster the remaining "
		<< "\n\thits.size() : " << hits.size()
		<< "\n\tinput\n\n";
    }

    
    auto michel_start = michel._start;
    
    for (const auto& thishit: hits)
      if (michel_start.SqDist(thishit) <= radius * radius)
	michel.push_back(thishit);
    
    

    michel._charge = 0;
    
    for(const auto& michel_hit : michel)
      michel._charge += michel_hit._q;
    
  }



}
#endif
