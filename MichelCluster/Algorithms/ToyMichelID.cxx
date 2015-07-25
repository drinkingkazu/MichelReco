#ifndef MICHELCLUSTER_TOYMICHELID_CXX
#define MICHELCLUSTER_TOYMICHELID_CXX

#include "ToyMichelID.h"

namespace michel {

  void ToyMichelID::EventReset()
  {}
  
  Michel ToyMichelID::Identify(const MichelCluster& cluster)
  { return Michel(); }

}
#endif
