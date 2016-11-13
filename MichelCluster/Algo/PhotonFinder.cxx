
#ifndef PHOTONFINDER_CXX
#define PHOTONFINDER_CXX

#include "PhotonFinder.h"
#include <sstream>
#include <cmath>
#include "Fmwk/ClusterVectorCalculator.h"

namespace michel {

  PhotonFinder::PhotonFinder()
  {
    _name = "PhotonFinder";
    _max_radius = 100.;
    _min_dot = 0.7;
  }

  void PhotonFinder::EventReset()
  {}
  
  bool PhotonFinder::ProcessCluster(MichelCluster& cluster,
				    const std::vector<HitPt>& hits){

    if (hits.size() == 0) return false;
    
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

    // get 2D direction of michel
    auto michelDir = cluster._michel._dir;

    // make sure direction vector is unit-normalized
    // necessary for subsequent calculations
    double mag = michelDir._w * michelDir._w + michelDir._t * michelDir._t;
    //std::cout << "Michel direction vector magnitude is " << mag << std::endl;

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
      ss << "michel start : [" << start._w << ", " << start._t << "]"
	 << "michel end   : [" << end._w   << ", " << end._t   << "]"
	 << "max radius   : " << sqrt(dMax);
      Print(msg::kINFO,this->Name(),ss.str());
    }

    // loop over all hits in the event
    // if within a reasonable radius [ 1 meter ]
    // then look further for the hit to be in a cone
    // starting from the michel start point
    // and pointing in the direction of the Michel

    size_t ctr =0;
    
    for (auto& h : hits){

      if (h._pl != 2)
	continue;

      // if hit already in Michel or Muon cluster -> ignore
      bool there = false;
      for (auto& mh : michel) { if(mh._id == h._id) { there = true; break;} }
      for (auto& ch : cluster._hits) { if(ch._id == h._id) { there = true; break;} }
      if (there)
	continue;
      
      if (start.SqDist(h) > 100. * 100.)
	continue;

      // get vector starting at Michel start and pointing to this point
      double start2pt_w = h._w - start._w;
      double start2pt_t = h._t - start._t;

      // normalize vector
      double start2pt_m = sqrt( start2pt_w * start2pt_w + start2pt_t * start2pt_t ); 
      start2pt_w /= start2pt_m;
      start2pt_t /= start2pt_m;

      // calculate dot-product with vector in Michel direction
      double dot = start2pt_w * michelDir._w + start2pt_t * michelDir._t;

      if (dot < _min_dot)
	continue;
      
      michel.push_back(h);
      ctr += 1;
    }

    //std::cout << "Added " << ctr << " hits" << std::endl;
    
    return true;
  }


}
#endif
