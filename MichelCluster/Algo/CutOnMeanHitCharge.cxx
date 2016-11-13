#ifndef CUTONMEANHITCHARGE_CXX
#define CUTONMEANHITCHARGE_CXX

#include "CutOnMeanHitCharge.h"

namespace michel{

  CutOnMeanHitCharge::CutOnMeanHitCharge()
  {
    _name    = "CutOnMeanHitCharge";
    _max_qavg = 400.;
  }

  bool CutOnMeanHitCharge::ProcessCluster(MichelCluster& cluster,
					  const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;
    
    //  get number of Michel hits and total charge in them
    double Qtot = 0.;
    int num_michel_hits = cluster._michel.size();

    for (size_t i=0; i < cluster._michel.size(); i++)
      Qtot += cluster._michel[i]._q;

    if (Qtot/num_michel_hits > _max_qavg)
      return false;

    return true;
  }

}

#endif
