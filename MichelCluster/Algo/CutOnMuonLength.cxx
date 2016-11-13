#ifndef CUTONMUONLENGTH_CXX
#define CUTONMUONLENGTH_CXX

#include "CutOnMuonLength.h"

namespace michel{

  CutOnMuonLength::CutOnMuonLength()
  {
    _name    = "CutOnMuonLength";
    _min_muon_length = 10;
  }

  bool CutOnMuonLength::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    // get the hit vector
    auto const& hit_v = cluster._hits;
    // get the boundary point
    auto const& boundary = (int)cluster._boundary;
    // is the michel forward (from boundary to end) or back
    auto const& forward = cluster._forward;
    // michel's start point
    auto const& start = cluster._hits[boundary];

    // ***************************************************
    // the boundary defines where the michel starts in the
    // ordered hit list. If less then 0 then wtf
    if (boundary < 0){
      return false;
    }

    // if the muon is too short -> ignore
    double muon_len = 0;
    if (forward){
      auto& mu_s = hit_v[0];
      muon_len = (mu_s._w-start._w)*(mu_s._w-start._w) + (mu_s._t-start._t)*(mu_s._t-start._t);
    }
    else{
      auto& mu_s = hit_v[hit_v.size()-1];
      muon_len = (mu_s._w-start._w)*(mu_s._w-start._w) + (mu_s._t-start._t)*(mu_s._t-start._t);
    }

    // *************************************************
    // require the muon segment to have a minimum length
    if (muon_len < _min_muon_length*_min_muon_length)
      return false;

    // if muon length is above cut -> proceed!
    return true;
  }

}

#endif
