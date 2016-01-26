#ifndef MICHELCLUSTER_CALCTRUNCATED_CXX
#define MICHELCLUSTER_CALCTRUNCATED_CXX

#include "CalcTruncated.h"
#include "Fmwk/MichelException.h"
#include "ClusterVectorCalculator.h"
#include <cmath>
#include <sstream>

namespace michel {

  CalcTruncated::CalcTruncated()
      :
      _covariance_window(11),
      _n_window_size(15),
      _p_above(0.25),
      _window_cutoff(3),
      _edgefix(3)
  {
    _name="CalcTruncated";
    _t0 = _t1 = _t2 = _t3 = 0;
    _n0 = _n1 = _n2 = _n3 = 0;

    return;
  }
  
  void CalcTruncated::EventReset()
  {

  }
  
  bool CalcTruncated::ProcessCluster(MichelCluster& cluster,
				     const std::vector<HitPt>& hits)
  { 
    
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
    
    //do truncated mean
    _watch.Start();
    truncated_mean = _clusterCalc.calc_smooth_mean(cluster,
						   _n_window_size,
						   _window_cutoff,
						   _p_above);
    _t0 += _watch.RealTime();
    _n0 += 1;
    
    if(truncated_mean.size() < _edgefix) {
      std::stringstream ss;
      ss << std::endl
	 << "\tUnable to fix edges on truncated mean, edgefix size: " << _edgefix
	 << "\t and truncated_mean.size(): " << truncated_mean.size() << std::endl;
      Print(msg::kERROR,__FUNCTION__,ss.str());
      return false;
    }
    
    for(int i = 0 ; i < _edgefix; ++i) {
      truncated_mean.at(i) = truncated_mean[_edgefix];
      truncated_mean.at(truncated_mean.size() - i - 1) = truncated_mean[truncated_mean.size() - _edgefix];
    }
    
    //Directionality considerations
    int dir_window = _covariance_window;
    _watch.Start();
    covariance     = _clusterCalc.calc_covariance(cluster._hits,dir_window);
    _t1 += _watch.RealTime();
    _n1 += 1;
    _watch.Start();
    slope          = _clusterCalc.calc_slope     (cluster._hits,dir_window);
    _t2 += _watch.RealTime();
    _n2 += 1;
    
    // must be odd, currently has no setter,
    // sorry that this method has no info on it, ask vic
    int s = 3; 
    _watch.Start();
    truncated_dqds = _clusterCalc.calc_smooth_derive(cluster._s_v,truncated_mean,s);
    _t3 += _watch.RealTime();
    _n3 += 1;
        
    //Lets play with truncated mean shaving...
    if(_verbosity <= msg::kINFO) {
      std::stringstream ss;
      ss << "\n\t\tIn CalcTruncated" << std::endl
	 << "\tI have " << truncated_mean.size() << " truncated mean size" << std::endl
	 << "\twith   " << truncated_dqds.size() << " derivative points." << std::endl
	 << "\tMy incoming cluster has " << cluster._hits.size() << " hits in it...";
      Print(msg::kINFO,__FUNCTION__,ss.str());
    }
    
    // put the informataion inside cluster object
    std::swap(cluster._t_mean_v,truncated_mean);
    std::swap(cluster._t_dqds_v,truncated_dqds);
    std::swap(cluster._chi2_v  ,covariance);
    std::swap(cluster._dirs_v  ,slope);

    return true;
  }

  void CalcTruncated::Report(){

    std::cout << "T0/clus = " << 1e6*_t0/(double)_n0 << std::endl
	      << "T1/clus = " << 1e6*_t1/(double)_n1 << std::endl
	      << "T2/clus = " << 1e6*_t2/(double)_n2 << std::endl
	      << "T3/clus = " << 1e6*_t3/(double)_n3 << std::endl;

  }
  

}
#endif
