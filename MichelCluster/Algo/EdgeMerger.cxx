#ifndef MICHELCLUSTER_EDGEMERGER_CXX
#define MICHELCLUSTER_EDGEMERGER_CXX

#include "EdgeMerger.h"

namespace michel{
  
  void EdgeMerger::EventReset()
  {}

  bool EdgeMerger::Merge(const MichelCluster& a, const MichelCluster& b)
  {
    /// Do not currently have confiure function yet, just use "6 [cm]" dist as before..."
    if(_verbosity <=  msg::kDEBUG) {
      std::cout << "Cluster A start/end:" << std::endl
		<< a._start.Print().c_str() << std::endl
		<< a._end.Print().c_str() << std::endl;
      std::cout << "Cluster B start/end:" << std::endl
		<< b._start.Print().c_str() << std::endl
		<< b._end.Print().c_str() << std::endl
		<< std::endl;
    }
    
    /// True if near, false if not
    return Touching(a,b,_edge_dist);
  }

  double EdgeMerger::Priority(const MichelCluster& cluster)
  {
    if(cluster._hits.size() < 2) return -1;

    return (double)(cluster._hits.size());
  }

  
  bool EdgeMerger::Touching (const MichelCluster& lhs,
			     const MichelCluster& rhs,
			     const double min_dist) const
  {
    return  ( lhs._start.Dist(rhs._end)   < min_dist ||
	      lhs._end.Dist  (rhs._start) < min_dist ||
	      lhs._start.Dist(rhs._start) < min_dist ||
	      lhs._end.Dist  (rhs._end)   < min_dist );
  }
  
}

#endif
