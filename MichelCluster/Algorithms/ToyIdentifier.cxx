#ifndef MICHELCLUSTER_TOYIDENTIFIER_CXX
#define MICHELCLUSTER_TOYIDENTIFIER_CXX

#include "ToyIdentifier.h"

namespace michel {

  void ToyIdentifier::EventReset()
  {}
  
  Michel ToyIdentifier::Identify(const MichelCluster& cluster)
  { return Michel(); }

}
#endif
