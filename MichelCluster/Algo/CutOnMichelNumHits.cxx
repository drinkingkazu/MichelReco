#ifndef CUTONMICHELNUMHITS_CXX
#define CUTONMICHELNUMHITS_CXX

#include "CutOnMichelNumHits.h"
#include <sstream>

namespace michel{

  CutOnMichelNumHits::CutOnMichelNumHits()
  {
    _name    = "CutOnMichelNumHits";
    _min_num_hits = 0;
    _max_num_hits = 100;
  }

  bool CutOnMichelNumHits::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    // get number of hits in the michel
    int num_michel_hits = cluster._michel.size();

    if (num_michel_hits < _min_num_hits){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "Num of michel hits is less than min allowed [" << num_michel_hits << "]";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }

    if (num_michel_hits > _max_num_hits){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "Num of michel hits is greater than min allowed [" << num_michel_hits << "]";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }

    // if muon length is above cut -> proceed!
    return true;
  }

}

#endif
