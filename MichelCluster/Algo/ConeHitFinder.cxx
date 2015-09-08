#ifndef CONEHITFINDER_CXX
#define CONEHITFINDER_CXX

#include "ConeHitFinder.h"
#include <sstream>
#include <cmath>
#include "ClusterVectorCalculator.h"

namespace michel {

  ConeHitFinder::ConeHitFinder()
  {
    _name = "ConeHitFinder";
    _max_radius = 20;
    _max_perp_distance = 5;
  }

  void ConeHitFinder::EventReset()
  {}
  
  bool ConeHitFinder::ProcessCluster(MichelCluster& cluster,
					   const std::vector<HitPt>& hits){

    // calculator tool
    ClusterVectorCalculator _clusterCalc;

    double dMax = _max_radius*_max_radius;

    auto& michel = cluster._michel;

    // How does this algorithm work:
    // calculate a linear fit to the michel points
    // for hits in the neighboring region
    // check if they are alligned with the slope
    // if so add them to the cluster

    
    //This michel was bogus when it came in, don't cluster further
    if (michel.size() == 0) return false;


    // fit michel to a straight line
    auto fitinfo = _clusterCalc.GetLinearFit(michel);

    if ( isnan(fitinfo.first) == true)
      return true;
    if (isnan(fitinfo.second) == true)
      return true;

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

    // get slope and intercept for linear fit to hits of michel cluster
    double slope = fitinfo.first;
    double intercept = fitinfo.second;

    // get end point as point projected forward
    // on line from michel start for a length
    // equal to the michel length
    double michelLen = sqrt(start.SqDist(end));
    double endW = start._w + michelLen;
    double endT = start._t + slope * michelLen;
    

    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "michel slope     : " << slope << "\t"
	 << "michel intercept : " << intercept << std::endl;
      Print(msg::kINFO,this->Name(),ss.str());
    }
    
    // loop over points and find:
    // 0) in front of in back of michel? if in back -> reject
    // 1) perpendicular distance to fit line

    size_t cnts = 0;
    for (auto h : nearbyHits){

      // if the point is backwards w.r.t. michel direction
      // ignore
      if ( ( (h._w-start._w)*(endW-start._w) +
	     (h._t-start._t)*(endT-start._t) ) < 0 )
	continue;

      // perpendicular distance of point to line:
      double d_perp = _clusterCalc.GetPerpendicularDistance(h,slope,intercept);
      double d_start = h.SqDist(start);

      if (d_perp < _max_perp_distance){
	michel.push_back(h);
	cnts += 1;
      }
    }
    
    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "hits added : " << cnts << std::endl;
      Print(msg::kINFO,this->Name(),ss.str());
    }

    return true;
  }



}
#endif
