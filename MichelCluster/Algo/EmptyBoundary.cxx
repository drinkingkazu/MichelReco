#ifndef MICHELCLUSTER_EMPTYBOUNDARY_CXX
#define MICHELCLUSTER_EMPTYBOUNDARY_CXX

#include "EmptyBoundary.h"

namespace michel {

void EmptyBoundary::EventReset()
{}

HitIdx_t EmptyBoundary::Boundary(MichelCluster& cluster) {

  return 0;

}

}
#endif
