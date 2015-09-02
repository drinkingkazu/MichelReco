#ifndef MICHELCLUSTER_EMPTYBOUNDARY_CXX
#define MICHELCLUSTER_EMPTYBOUNDARY_CXX

#include "EmptyBoundary.h"

namespace michel {

  void EmptyBoundary::EventReset(){}
  
  void EmptyBoundary::ProcessCluster(MichelCluster& cluster,
				     const std::vector<HitPt>& hits)
  {
    return;
  }
  
}
#endif
