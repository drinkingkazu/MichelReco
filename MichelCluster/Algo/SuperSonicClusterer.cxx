#ifndef SUPERSONICCLUSTERER_CXX
#define SUPERSONICCLUSTERER_CXX

#include "SuperSonicClusterer.h"
#include <sstream>
#include <cmath>

namespace michel {

  SuperSonicClusterer::SuperSonicClusterer()
  {
    _merge_till_converge = false;
    _name = "SuperSonicClusterer";
    _max_radius = 10;
    _hit_radius = 3;
    _use_hit_radius = false;
  }

  void SuperSonicClusterer::EventReset()
  {}
  
  bool SuperSonicClusterer::ProcessCluster(MichelCluster& cluster,
					   const std::vector<HitPt>& hits){

    if (hits.size() == 0) return false;

    double dMax = _max_radius*_max_radius;

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
    double lenMax = 0;
    for (auto& h : michel){
      double dist = start.SqDist(h);
      if (dist > lenMax){
	lenMax = dist;
	end = h;
      }// if the distance is largest
    }// for all hits in michel
	
    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "michel start : [" << start._w << ", " << start._t << "]" << std::endl
	 << "michel end   : [" << end._w   << ", " << end._t   << "]"  << std::endl
	 << "max radius   : " << sqrt(dMax);
      Print(msg::kINFO,this->Name(),ss.str());
    }


    // get a list of hits that is less than dMax away from either the start or end point
    std::vector<michel::HitPt> nearbyHits;
    for (auto& h : hits){

      if (h._pl != 2) continue;

      bool there = false;
      for (auto& mh : michel) { if(mh._id == h._id) { there = true; break;} }
      for (auto& ch : cluster._hits) { if(ch._id == h._id) { there = true; break;} }

      if (there) continue;

      if ( (start.SqDist(h) < dMax) or (end.SqDist(h) < dMax) )
	nearbyHits.push_back(h);
    }

    // merge?
    bool merge = true;
    int step = 0;
    while(merge){

      if (_verbosity <= msg::kINFO){
	std::stringstream ss;
	ss << "merging step " << step << "\t nearbyHits size : " << nearbyHits.size();
	Print(msg::kINFO,this->Name(),ss.str());
      }
      
      
      // get a list of hits that is less than dMax away from either the start or end point
      std::vector<michel::HitPt> nearbyHits_cpy;
      bool there;
      for (auto& h : nearbyHits) {
	
	if (h._pl != 2) continue;
	
	there = false;
	
	for (auto& mh : michel) { if(mh._id == h._id) { there = true; break;} }
	if(there) continue;
	for (auto& ch : cluster._hits) { if(ch._id == h._id) { there = true; break;} }
	if(there) continue;
	
	if ( (start.SqDist(h) < dMax) or (end.SqDist(h) < dMax) )
	  nearbyHits_cpy.push_back(h);
      }

      nearbyHits.clear();
      nearbyHits = nearbyHits_cpy;

      
      // now, for each hit in the michel
      // 1) determine its distance from the start point
      // 2) loop through all hits in the "nearbyHits" vector
      // 3) if any of them are closer than the start point is to
      // the specific michel hit we are at, include them
      // have a vector where additional hits are temporarily stored
      std::vector<michel::HitPt> newHits;
      for (auto& h : nearbyHits){
	for (auto& mh : michel){
	  // if we should use a fixed hit radius
	  if (_use_hit_radius){
	    if (mh.SqDist(h) < _hit_radius){
	      newHits.push_back(h);
	      break;
	    }
	  }
	  // else, if we use the distance to the start
	  else{
	    // distance to start:
	    auto dstart = start.SqDist(mh);
	    if (mh.SqDist(h) < dstart){
	      newHits.push_back(h);
	      break;
	    }
	  }// if we should use distance to start
	}// for all michel hits
      }// for all nearby hits

      if (_verbosity <= msg::kINFO){
	std::stringstream ss;
	ss << "added " << newHits.size() << " hits";
      Print(msg::kINFO,this->Name(),ss.str());
      }
      
      // now that we found potential new hits, append them to the michel object
      for (auto &h : newHits)
	michel.push_back(h);
      
      
      // recalculte the michel's charge taking into account the new hits
      michel._charge = 0;
      for(const auto& michel_hit : michel)
	michel._charge += michel_hit._q;
    
      // if we should merge-till-converge -> keep merge true
      // if we have not added any new hits -> turn merge off
      if ( (newHits.size() == 0) or (!_merge_till_converge) )
	merge = false;
      
      step += 1;
    }// while merge is true
    
    return true;
  }



}
#endif
