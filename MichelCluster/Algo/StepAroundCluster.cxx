#ifndef STEPAROUNDCLUSTER_CXX
#define STEPAROUNDCLUSTER_CXX

#include "StepAroundCluster.h"
#include <algorithm>

namespace michel {

  void StepAroundCluster::EventReset()
  {}
  
  bool StepAroundCluster::ProcessCluster(MichelCluster& cluster,
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
    }
    
    michel._charge = 0;
    
    for(const auto& michel_hit : michel)
      michel._charge += michel_hit._q;

    return true;
  }
  
  

}
#endif
