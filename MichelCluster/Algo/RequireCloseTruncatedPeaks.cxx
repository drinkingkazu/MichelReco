#ifndef MICHELCLUSTER_REQUIRECLOSETRUNCATEDPEAKS_CXX
#define MICHELCLUSTER_REQUIRECLOSETRUNCATEDPEAKS_CXX

#include "RequireCloseTruncatedPeaks.h"
#include "Fmwk/MichelException.h"
#include "Fmwk/ClusterVectorCalculator.h"
#include <cmath>
#include <sstream>

namespace michel {
  
  void RequireCloseTruncatedPeaks::EventReset()
  {}
  
  bool RequireCloseTruncatedPeaks::ProcessCluster(MichelCluster& cluster,
						  const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;
    
    if(cluster._t_mean_v.size() < cluster._hits.size()) {
      Print(msg::kEXCEPTION,this->Name(),"Truncated Mean vector size less than num hits, run CalcTruncated");
      throw MichelException();
    }
    if(cluster._t_dqds_v.size() < cluster._hits.size()) {
      Print(msg::kEXCEPTION,this->Name(),"dQdS vector size less than num hits, run CalcTruncated");
      throw MichelException();
    }

    ClusterVectorCalculator _clusterCalc;
    
    //With this new information, calculate the boundary point between possible muon end and michel start
    
    size_t candidate_loc     = _clusterCalc.find_max(cluster._t_mean_v);
    size_t dqdscandidate_loc = _clusterCalc.find_min(cluster._t_dqds_v); 
    
    if((candidate_loc     >= cluster._hits.size()))
      return false;
    
    if((dqdscandidate_loc >= cluster._hits.size()))
      return false;
    
    if(abs(dqdscandidate_loc - candidate_loc) > _maxDistance)
      return false;

    return true;
  }

}
#endif
