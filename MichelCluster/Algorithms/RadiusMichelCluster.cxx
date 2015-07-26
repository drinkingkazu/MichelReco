#ifndef RADIUSMICHELCLUSTER_CXX
#define RADIUSMICHELCLUSTER_CXX

#include "RadiusMichelCluster.h"

namespace michel {

  void RadiusMichelCluster::EventReset()
  {}
  
  void RadiusMichelCluster::Cluster(Michel& michel,
				    const std::vector<HitPt>& hits){

    if (michel.size() == 0) return;

    double min_radius = 10; //hardcoded for now
    double radius = michel._length;
    
    //set to min rad
    if (radius < min_radius)
      radius = min_radius;
    
    
    
    HitPt michel_start = michel._start; //probably want to move to end...

    for (const auto& thishit: hits){
      //if within circle, include in michel
      if (michel_start.SqDist(thishit) <= radius * radius){
	michel.push_back(thishit);
      }
    }
  }


}
#endif
