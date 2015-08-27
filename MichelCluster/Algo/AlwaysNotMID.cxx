#ifndef ALWAYSNOTMID_CXX
#define ALWAYSNOTMID_CXX

#include "AlwaysNotMID.h"

namespace michel{

  AlwaysNotMID::AlwaysNotMID() { }

  bool AlwaysNotMID::IsMichel(const MichelCluster& michel,
				      const std::vector<HitPt>& hits)
  {
    return true; /// It's a michel always
  }
  
}

#endif
