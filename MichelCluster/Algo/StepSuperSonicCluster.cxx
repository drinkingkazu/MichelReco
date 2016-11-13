#ifndef STEPSUPERSONICCLUSTER_CXX
#define STEPSUPERSONICCLUSTER_CXX

#include "StepSuperSonicCluster.h"
#include <algorithm>
#include <sstream>
#include <cmath>

namespace michel {

  StepSuperSonicCluster::StepSuperSonicCluster()
  {
    _name = "StepSuperSonicCluster";
    _merge_till_converge = false;
    _max_radius = 10;
    _hit_radius = 3;
    _use_hit_radius = false;
  }

  void StepSuperSonicCluster::EventReset()
  {}
  
  bool StepSuperSonicCluster::ProcessCluster(MichelCluster& cluster,
					     const std::vector<HitPt>& hits){

    if (hits.size() == 0) return false;
    
    auto& michel = cluster._michel;
    
    /// This michel was bogus when it came in, don't cluster further
    if (michel.size() == 0) return false;
    
    auto michel_end  = michel[michel.size() - 1];
    auto max_step_sq = _max_step * _max_step;
    
    std::vector<size_t> taken_hits;
    taken_hits.reserve(hits.size());

    
    /// Step around on the end point a bit, jumping to the nearest
    /// hit with cutoff < _max_step.


    //temporary fix for SSS is add cluster hit ids to taken vector
    for (auto& ch : cluster._hits) taken_hits.push_back(ch._id);

    if (_verbosity == msg::kINFO)
      std::cout << "michel size before step: " << michel.size() << std::endl;
    
    while(1) {
      
      bool   done = true;
      double dist = kINVALID_DOUBLE;
      HitPt  close;
      
      for (const auto& thishit: hits) {

	/// if I already added this hit move on, std::find probably killer slow, oh well
	if(std::find(taken_hits.begin(),
		     taken_hits.end(),
		     thishit._id) != taken_hits.end()) continue;
	
	auto curr_dist = michel_end.SqDist(thishit);
	if ( curr_dist < dist && curr_dist < max_step_sq) {
	  close = thishit;
	  dist = curr_dist;
	  done = false;
	}
      }

      if(!done) {
	michel    .push_back(close);
	taken_hits.push_back(close._id);
       
	michel_end  = michel[michel.size() - 1];
      }
      
      if(done) break;

      //return true;
    }

    if (_verbosity == msg::kINFO)
      std::cout << "michel size after step: " << michel.size() << std::endl;
    
    //return true;
    
    
    /// Add on David C algo, after we step around and make something, we then super sonic based
    /// on distance between michel start and location in michel, any other hits we find
    /// in this vicinity we to michel in no particular order
    
    /// some input hits are now both michel and input, I don't know what's better at this point
    /// if we don't care about speed then I can void out the ones that were taken before
    /// or maybe since the michel is small might as well loop over it and check ID to avoid using
    /// std::find over and over

    double dMax = _max_radius*_max_radius;

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
	ss << "merging step " << step;
	Print(msg::kINFO,this->Name(),ss.str());
      }
      
      
      // get a list of hits that is less than dMax away from either the start or end point
      std::vector<michel::HitPt> nearbyHits_cpy;
      bool there;
      for (auto& h : nearbyHits) {
	
	if (h._pl != 2) continue;
	
	there = false;
	
	for (auto& mh : michel) { if(mh._id == h._id) { there = true; break;} }
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
      for (auto& mh : michel){
	// if we should use a fixed hit radius
	if (_use_hit_radius){
	  for (auto& h : nearbyHits){
	    if (mh.SqDist(h) < _hit_radius)
	      newHits.push_back(h);
	  }
	}
	// else, if we use the distance to the start
	else{
	  // distance to start:
	  auto dstart = start.SqDist(mh);
	  for (auto& h : nearbyHits){
	    if (mh.SqDist(h) < dstart)
	      newHits.push_back(h);
	  }
	}// if we should use distance to start
      }// for all michel hits
      
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
      if (_verbosity <= msg::kINFO){
	std::stringstream ss;
	ss << "added " << newHits.size() << " hits";
	Print(msg::kINFO,this->Name(),ss.str());
      }
      if ( (newHits.size() == 0) or (!_merge_till_converge) )
	merge = false;
      
      step += 1;
    }// while merge is true
    
    return true;
  }

    

}
#endif
