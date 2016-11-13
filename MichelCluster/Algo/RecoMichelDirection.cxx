#ifndef RECOMICHELDIRECTION_CXX
#define RECOMICHELDIRECTION_CXX

#include "RecoMichelDirection.h"
#include <sstream>

namespace michel{

  RecoMichelDirection::RecoMichelDirection()
  {
    _name    = "RecoMichelDirection";
  }

  bool RecoMichelDirection::ProcessCluster(MichelCluster& cluster,
					  const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;
    
    // grab hits from hits associated to Michel electron part only
    // (no photons) and find a 2D direction for this Michel

    auto const& michel_hits          = cluster._michel;
    auto const& electron_hit_indices = cluster._michel._electron_hit_idx_v;
    auto const& michel_start         = cluster._michel._start;

    // are there Michel hits? if no -> quit
    if (electron_hit_indices.size() == 0){
      if(_verbosity >= msg::kWARNING) {
	std::stringstream ss;
	ss << "\n\t\t NO MICHEL HITS !?!? " << std::endl;
	Print(msg::kERROR,__FUNCTION__,ss.str());
      }
      return false;
    }

    // get vector of electron hits
    std::vector<michel::HitPt> elecPts;
    
    for (auto const& idx : electron_hit_indices)
      elecPts.push_back( michel_hits[idx] );
    
    // calculate 2D direction of this cluster [charge weighted]
    HitPt dir(0,0,0,0,0);

    double weightedDir_w = 0;
    double weightedDir_t = 0;
    double Qtot = 0;

    for (auto const& hit : elecPts){

      double dt = hit._t - michel_start._t;
      double dw = hit._w - michel_start._w;

      weightedDir_t += dt * hit._q;
      weightedDir_w += dw * hit._q;
      Qtot += hit._q;

    }

    weightedDir_t /= Qtot;
    weightedDir_w /= Qtot;

    // normalize length
    double mag = sqrt( weightedDir_w * weightedDir_w + weightedDir_t * weightedDir_t );

    cluster._michel._dir = HitPt( Qtot, weightedDir_w / mag, weightedDir_t / mag, 0, 2);

    return true;
  }

}

#endif
