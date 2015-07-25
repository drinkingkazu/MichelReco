#ifndef TOYMICHELCLUSTER_CXX
#define TOYMICHELCLUSTER_CXX

#include "ToyMichelCluster.h"

namespace michel {

  void ToyMichelCluster::EventReset()
  {}
  
  Michel ToyMichelCluster::Cluster(const MichelCluster& cluster,
				   const std::vector<HitPt>& hits)
  { return cluster._michel; }

}
#endif
