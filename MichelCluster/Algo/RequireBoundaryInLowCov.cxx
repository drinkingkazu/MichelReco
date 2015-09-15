#ifndef MICHELCLUSTER_REQUIREBOUNDARINLOWCOV_CXX
#define MICHELCLUSTER_REQUIREBOUNDARINLOWCOV_CXX

#include "RequireBoundaryInLowCov.h"
#include "Fmwk/MichelException.h"
#include <cmath>
#include <sstream>
namespace michel {
  
  void RequireBoundaryInLowCov::EventReset()
  {}
  
  bool RequireBoundaryInLowCov::ProcessCluster(MichelCluster& cluster,
					       const std::vector<HitPt>& hits)
  { 

    if(cluster._chi2_v.size() < cluster._hits.size()) {
      Print(msg::kEXCEPTION,this->Name(),"Covariance vector size less than num hits, run CalcTruncated");
      throw MichelException();
    }
    
    const auto& covariance = cluster._chi2_v;
    const auto& idx = cluster._boundary;

    // move on only if avg. covariance for hits near the "start" is less than some value
    // idx is the position of the reconstructed michel start
    double avg_covariance = 0;
    int    counts = 0;

    for (int i = (idx-3); i < (idx+4); i++){ // -3 and +4 is hardcoded ok fine
      // make sure we fall in the vector's range
      if ( (i >= 0) and (i < covariance.size()) ){
	avg_covariance += fabs(covariance[i]);
	counts += 1;
      }
    }

    // make sure we have at least 1 point!
    if (avg_covariance == 0)
      return false;

    // if so, check that the average covariance is below
    // the max allowed covariance
    avg_covariance /= counts;

    if (avg_covariance > _maxCovarianceAtStart)
      return false;

    return true;
  }

}
#endif
