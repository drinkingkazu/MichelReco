#ifndef MICHELCLUSTER_TRUNCATEDQBOUNDARY_CXX
#define MICHELCLUSTER_TRUNCATEDQBOUNDARY_CXX

#include "TruncatedQBoundary.h"
#include "Fmwk/MichelException.h"
#include "Fmwk/ClusterVectorCalculator.h"
#include <cmath>
#include <sstream>
namespace michel {
  
  void TruncatedQBoundary::EventReset()
  {}
  
  bool TruncatedQBoundary::ProcessCluster(MichelCluster& cluster,
					  const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    /// call instance of "ClusterVectorCalculator"
    /// this class has a bunch of utility functions
    /// to calculate stuff based on the vecor of
    /// hits for a cluster
    ClusterVectorCalculator _clusterCalc;
    
    std::vector<double> truncated_mean;
    std::vector<double> truncated_dqds;
    std::vector<double> covariance; 
    std::vector<double> slope; 
    
    truncated_mean.reserve(cluster._ordered_pts.size());
    truncated_dqds.reserve(cluster._ordered_pts.size());
    covariance.reserve    (cluster._ordered_pts.size());
    slope.reserve         (cluster._ordered_pts.size());
    
    //hardcoded for now will become configurable
    // double _n_window_size = 15;
    // double _p_above       = 0.25;
    // int    _window_cutoff = 3;

    //do truncated mean
    truncated_mean = _clusterCalc.calc_smooth_mean(cluster,
						   _n_window_size,
						   _window_cutoff,
						   _p_above);

    if(truncated_mean.size() < _edgefix) {
      std::stringstream ss;
      ss << "\tUnable to fix edges on truncated mean, edgefix size: " << _edgefix
	 << "\t and truncated_mean.size(): " << truncated_mean.size();
      Print(msg::kERROR,__FUNCTION__,ss.str());
      return false;
    }
    
    for(int i = 0 ; i < _edgefix; ++i) {
      truncated_mean.at(i) = truncated_mean[_edgefix];
      truncated_mean.at(truncated_mean.size() - i - 1) = truncated_mean[truncated_mean.size() - _edgefix];
    }
    
    //Directionality considerations
    int dir_window = _covariance_window;
    covariance     = _clusterCalc.calc_covariance(cluster._hits,dir_window);
    slope          = _clusterCalc.calc_slope     (cluster._hits,dir_window);
    
    // must be odd, currently has no setter,
    // sorry that this method has no info on it, ask vic
    int s = 3; 
    truncated_dqds = _clusterCalc.calc_smooth_derive(cluster._s_v,truncated_mean,s);
        
    //Lets play with truncated mean shaving...
    if(_verbosity <= msg::kINFO) {
      std::stringstream ss;
      ss << "\t\tIn TruncatedQBoundary" << std::endl
	 << "\tI have " << truncated_mean.size() << " truncated mean size" << std::endl
	 << "\twith   " << truncated_dqds.size() << " derivative points." << std::endl
	 << "\tMy incoming cluster has " << cluster._hits.size() << " hits in it...";
      Print(msg::kINFO,__FUNCTION__,ss.str());
    }
    
    //With this new information, calculate the boundary point between possible muon end and michel start
    
    int candidate_loc     = _clusterCalc.find_max(truncated_mean);
    int dqdscandidate_loc = _clusterCalc.find_min(truncated_dqds); 

    std::swap(cluster._t_mean_v,truncated_mean);
    std::swap(cluster._t_dqds_v,truncated_dqds);
    
    if((candidate_loc     >= cluster._hits.size()))
      return false;
    
    if((dqdscandidate_loc >= cluster._hits.size()))
      return false;
    
    if(abs(dqdscandidate_loc - candidate_loc) > _maxDistance)
      return false;

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
    
    std::swap(cluster._chi2_v,covariance);
    std::swap(cluster._dirs_v,slope);

    cluster._boundary = cluster._ordered_pts[idx];

    return true;
  }
  
}
#endif
