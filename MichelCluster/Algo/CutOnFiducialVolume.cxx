#ifndef MICHEL_CUTONFIDUCIALVOLUME_CXX
#define MICHEL_CUTONFIDUCIALVOLUME_CXX

#include "CutOnFiducialVolume.h"
#include <sstream>

namespace michel {

CutOnFiducialVolume::CutOnFiducialVolume()
{
  _name    = "CutOnFiducialVolume";

}

bool CutOnFiducialVolume::ProcessCluster(MichelCluster& cluster,
    const std::vector<HitPt>& hits)
{

  // get number of hits in the michel
  int num_michel_hits = cluster._michel.size();

  if (!num_michel_hits) return false;

  // if muon length is above cut -> proceed!
  return true;
}

}

#endif
