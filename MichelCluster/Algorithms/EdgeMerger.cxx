#ifndef MICHELCLUSTER_EDGEMERGER_CXX
#define MICHELCLUSTER_EDGEMERGER_CXX

#include "EdgeMerger.h"

namespace michel{

  void EdgeMerger::EventReset()
  {}

  bool EdgeMerger::Merge(const MichelCluster& a, const MichelCluster& b)
  {
    
    return true;
  }

  double EdgeMerger::Priority(const MichelCluster& cluster)
  {
    if(cluster._hits.size() < 2) return -1;

    return (double)(cluster._hits.size());
  }
  
}

#endif
