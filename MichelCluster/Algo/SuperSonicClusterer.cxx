#ifndef SUPERSONICCLUSTERER_CXX
#define SUPERSONICCLUSTERER_CXX

#include "SuperSonicClusterer.h"

namespace michel {

  void SuperSonicClusterer::EventReset()
  {}
  
  bool SuperSonicClusterer::ProcessCluster(MichelCluster& cluster,
					   const std::vector<HitPt>& hits){

    auto& michel = cluster._michel;

    // How does this algorithm work:
    // for each hit along the michel,
    // we search in a circle of radius = distance from that hit to the michel's start point
    // any hit within that radius becomes part of the michel cluster
    // you can visualize this as searching in an area that has
    // a cone-like shape starting at the michel start point
    // hits close to the start point will have a small radius
    // hits further away will have a larger radius to search in

    
    //This michel was bogus when it came in, don't cluster further
    if (michel.size() == 0) return false;


    // first get the maximum distance:
    // point furthest from the start point defines
    // the maximum radius
    // get a hit-list of hits within that distance
    // of the start-point or end-point

    auto start = michel._start;
    // search for an end-point by looking for
    // the hit furthest from the michel start
    HitPt end;
    double dMax = 0;
    for (auto& h : michel){
      double dist = start.SqDist(h);
      if (dist > dMax){
	dMax = dist;
	end = h;
      }// if the distance is largest
    }// for all hits in michel
	
    if (_verbosity == msg::kINFO){
      std::cout << "michel start : [" << start._w << ", " << start._t << "]" << std::endl
		<< "michel end   : [" << end._w   << ", " << end._t   << "]"  << std::endl
		<< "max radius   : " << dMax << std::endl;
    }

    // get a list of hits that is less than dMax away from either the start or end point
    std::vector<michel::HitPt> nearbyHits;
    for (auto& h : hits){
      if ( (start.SqDist(h) < dMax) or (end.SqDist(h) < dMax) )
	nearbyHits.push_back(h);
    }

    // now, for each hit in the michel
    // 1) determine its distance from the start point
    // 2) loop through all hits in the "nearbyHits" vector
    // 3) if any of them are closer than the start point is to
    // the specific michel hit we are at, include them
    // have a vector where additional hits are temporarily stored
    std::vector<michel::HitPt> newHits;
    for (auto& mh : michel){
      // distance to start:
      auto dstart = start.SqDist(mh);
      for (auto& h : nearbyHits){
	if (mh.SqDist(h) < dstart)
	  newHits.push_back(h);
      }// for all nearby hits
    }// for all michel hits

    // now that we found potential new hits, append them to the michel object
    for (auto &h : newHits)
      michel.push_back(h);
    
    
    // recalculte the michel's charge taking into account the new hits
    michel._charge = 0;
    for(const auto& michel_hit : michel)
      michel._charge += michel_hit._q;
   

    return true;
  }



}
#endif
