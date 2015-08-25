#ifndef MICHELCLUSTER_MATCHQBOUNDARIES_CXX
#define MICHELCLUSTER_MATCHQBOUNDARIES_CXX

#include "MatchBoundaries.h"
#include "Fmwk/MichelException.h"
#include "ClusterVectorCalculator.h"
#include <cmath>

namespace michel {
  
  void MatchBoundaries::EventReset()
  {}
  
  HitIdx_t MatchBoundaries::Boundary(MichelCluster& cluster)
  { 

    /// call instance of "ClusterVectorCalculator"
    /// this class has a bunch of utility functions
    /// to calculate stuff based on the vecor of
    /// hits for a cluster
    ClusterVectorCalculator _clusterCalc;
    
    std::vector<double> truncated_mean;
    std::vector<double> truncated_dqds;
    std::vector<double> covariance; 
    
    truncated_mean.reserve(cluster._ordered_pts.size());
    truncated_dqds.reserve(cluster._ordered_pts.size());
    covariance.reserve    (cluster._ordered_pts.size());
    
    //hardcoded for now will become configurable
    double _n_window_size = 15;
    double _p_above       = 0.25;
    int    _window_cutoff = 3;

    //do truncated mean
    truncated_mean = _clusterCalc.calc_smooth_mean(cluster,
						   _n_window_size,
						   _window_cutoff,
						   _p_above);
    
    for(int i = 0 ; i < 3; ++i) {
      truncated_mean.at(i) = truncated_mean[3];
      truncated_mean.at(truncated_mean.size() - i - 1) = truncated_mean[truncated_mean.size() - 3];
    }
    
    //Directionality considerations
    int dir_window = 11;
    covariance     = _clusterCalc.calc_covariance(cluster._hits,dir_window);
    
    // must be odd, currently has no setter,
    // sorry that this method has no info on it, ask vic
    int s = 3; 
    truncated_dqds = _clusterCalc.calc_smooth_derive(cluster._s_v,truncated_mean,s);
        
    //Lets play with truncated mean shaving...
    if(_verbosity <= msg::kINFO) {
      std::cout << "\n\t\tIn MatchBoundaries\n"
		<< "\tI have " << truncated_mean.size() << " truncated mean size\n"
		<< "\twith   " << truncated_dqds.size() << " derivative points.\n"
		<< "\tMy incoming cluster has " << cluster._hits.size() << " hits in it...\n";
    }
    
    //With this new information, calculate the boundary point between possible muon end and michel start
    
    auto candidate_loc     = _clusterCalc.find_max(truncated_mean);
    auto dqdscandidate_loc = _clusterCalc.find_min(truncated_dqds); 

    std::swap(cluster._t_mean_v,truncated_mean);
    std::swap(cluster._t_dqds_v,truncated_dqds);
    
    if((candidate_loc     >= cluster._hits.size()))
      return kINVALID_SIZE;
    
    if((dqdscandidate_loc >= cluster._hits.size()))
      return kINVALID_SIZE;
    
    if(abs(dqdscandidate_loc - candidate_loc) > _maxDistance)
      return kINVALID_SIZE;
    

    auto right = cluster._ordered_pts.size() - 1 - candidate_loc;
    auto left  = candidate_loc;
    
    int  iMin = 0;
    int  iMax = 0;
    
    if(right >= _maxDistance) iMax  = _maxDistance   + candidate_loc;
    if(left  >= _maxDistance) iMin  = candidate_loc - _maxDistance;

    if(right < _maxDistance)  iMax  = cluster._hits.size() - 1;
    if(left  < _maxDistance)  iMin  = 0;

    // holder for hit with largest charge -> this will identify the location
    // of the hit that defines the michel boundary
    auto k   = 0.0;
    auto idx = 0;
    
    for(int w = iMin; w <= iMax; ++w) {
      auto c = cluster._hits[cluster._ordered_pts[w]]._q;
      // if this hit has more charge than any other
      if(c > k) { k = c; idx = w; }
    }

    // move on only if avg. covariance for hits near the "start" is less than some value
    // idx is the position of the reconstructed michel start
    double avg_covariance = 0;
    int counts = 0;
    for (int i = (idx-3); i < (idx+4); i++){
      // make sure we fall in the vector's range
      if ( (i >= 0) and (i < covariance.size()) ){
	avg_covariance += fabs(covariance[i]);
	counts += 1;
      }
    }
    // make sure we have at least 1 point!
    if (avg_covariance == 0)
      return kINVALID_SIZE;
    // if so, check that the average covariance is below
    // the max allowed covariance
    avg_covariance /= counts;
    if (avg_covariance > _maxCovarianceAtStart)
      return kINVALID_SIZE;
    
    std::swap(cluster._chi2_v,covariance);

    return cluster._ordered_pts[idx];
  }
  
}
#endif
