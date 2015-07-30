#ifndef MICHELCLUSTER_TOYBOUNDARY_CXX
#define MICHELCLUSTER_TOYBOUNDARY_CXX

#include "ToyBoundary.h"

namespace michel {

  void ToyBoundary::EventReset()
  {}
    
  HitIdx_t ToyBoundary::Boundary(MichelCluster& cluster)
  { return 0; }

}
#endif
