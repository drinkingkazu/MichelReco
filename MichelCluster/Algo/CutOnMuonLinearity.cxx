#ifndef CUTONMUONLINEARITY_CXX
#define CUTONMUONLINEARITY_CXX

#include "CutOnMuonLinearity.h"

namespace michel{

  CutOnMuonLinearity::CutOnMuonLinearity()
  {
    _name    = "CutOnMuonLinearity";
    _chi_min = 0.8;
    _frac_min_hits = 0.5;
  }

  bool CutOnMuonLinearity::ProcessCluster(MichelCluster& cluster,
					  const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    // get the chi vector
    auto const& chi_v = cluster._chi2_v;
    // get the hit vector
    auto const& hit_v = cluster._hits;
    // get the boundary point
    auto const& boundary = (int)cluster._boundary;
    // is the michel forward (from boundary to end) or back
    auto const& forward = cluster._forward;

    // ***************************************************
    // the boundary defines where the michel starts in the
    // ordered hit list. If less then 0 then wtf
    if (boundary < 0){
      return false;
    }

    // First start with some initial cuts
    // these are (stringent) sanity checks
    // to remove generally very bad michels

    // chi values for entire muon track
    std::vector<double> chi_muon;

    // if the muon is too short -> ignore
    if (forward){
      for (int i = 0; i < boundary; i++)
	chi_muon.push_back(fabs(chi_v[i]));
    }
    else{
      for (size_t i = boundary; i < hit_v.size()-1; i++)
	chi_muon.push_back(fabs(chi_v[i]));
    }

    // ******************************************************************
    // if < 50% of chi_v entries belonging to muon are above 80% -> remove
    int num_above_cut = 0;
    for (auto& c : chi_muon)
      if (c > _chi_min) num_above_cut += 1;
    double frac_above = ((double)num_above_cut) / ((double)chi_muon.size());
    if (frac_above < _frac_min_hits)
      return false;


    // if the muon is considered sufficiently linear -> return true!
    return true;
  }

}

#endif
