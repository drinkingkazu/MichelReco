#ifndef EMPTYMIDFILTER_CXX
#define EMPTYMIDFILTER_CXX

#include "EmptyMIDFilter.h"

namespace michel {

  EmptyMIDFilter::EmptyMIDFilter(){_name="EmptyMIDFilter";}
  
  bool EmptyMIDFilter::ProcessCluster(MichelCluster& cluster,
				      const std::vector<HitPt>& hits)
  {
    return true;
  }
  
}

#endif
