#ifndef CUTONTOTNUMHITS_CXX
#define CUTONTOTNUMHITS_CXX

#include "CutOnTotNumHits.h"
#include <sstream>

namespace michel{

  CutOnTotNumHits::CutOnTotNumHits()
  {
    _name    = "CutOnTotNumHits";
    _min_num_hits = 0;
  }

  bool CutOnTotNumHits::ProcessCluster(MichelCluster& cluster,
				       const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    // get number of hits in the cluster
    int num_hits = cluster._hits.size();

    if (num_hits < _min_num_hits){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "Num of cluster hits is less than min allowed [" << num_hits << "]";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }

    return true;
  }

}

#endif
