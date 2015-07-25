#ifndef MICHELCLUSTER_EDGEMERGER_CXX
#define MICHELCLUSTER_EDGEMERGER_CXX

#include "EdgeMerger.h"

namespace michel{

  void EdgeMerger::EventReset()
  {}

  bool EdgeMerger::Merge(const MichelCluster& a, const MichelCluster& b)
  {

    //Do not currently have confiure function yet, just use "6 [cm]" dist as before..."
    if(a.Touching(b,6)) { return true;  }
    else                { return false; }
  }

  double EdgeMerger::Priority(const MichelCluster& cluster)
  {
    if(cluster._hits.size() < 2) return -1;

    return (double)(cluster._hits.size());
  }
  
}

#endif
