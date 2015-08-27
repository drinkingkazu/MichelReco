#ifndef EMPTYMIDFILTER_CXX
#define EMPTYMIDFILTER_CXX

#include "EmptyMIDFilter.h"

namespace michel {

EmptyMIDFilter::EmptyMIDFilter()
{

}

bool EmptyMIDFilter::IsMichel(const MichelCluster& michel,
                              const std::vector<HitPt>& hits)
{
    return true;
}

}

#endif
