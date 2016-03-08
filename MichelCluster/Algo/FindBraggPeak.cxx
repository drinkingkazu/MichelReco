#ifndef FINDBRAGGPEAK_CXX
#define FINDBRAGGPEAK_CXX

#include "FindBraggPeak.h"
#include "Fmwk/ClusterVectorCalculator.h"
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>

namespace michel {

  FindBraggPeak::FindBraggPeak()
  {
    _name = "FindBraggPeak";
    _minBraggArea = 0.;
  }

  void FindBraggPeak::EventReset()
  {}
  
  bool FindBraggPeak::ProcessCluster(MichelCluster& cluster,
					   const std::vector<HitPt>& hits){

    /// call instance of "ClusterVectorCalculator"
    /// this class has a bunch of utility functions
    /// to calculate stuff based on the vecor of
    /// hits for a cluster
    ClusterVectorCalculator _calc;

    // Ths algorithm calculates the area of the Bragg Peak
    // for a cluster
    // 1) the bragg peak amplitude and dS position is found
    //    using the truncated mean vector
    std::pair<size_t,double> braggInfo = _calc.GetMaxIndex(cluster._t_mean_v);
    // set the bragg position for the cluster boundary
    HitIdx_t braggBoundary = braggInfo.first;
    // if we find the bragg peak in the 1st half
    // -> flip the order of the vector so that
    // the muon is moving towards the right
    if ( braggInfo.first < (size_t)(cluster._t_mean_v.size()/2.) ){
      std::reverse(cluster._s_v.begin(),cluster._s_v.end());
      std::reverse(cluster._t_mean_v.begin(),cluster._t_mean_v.end());
      size_t newIdx = cluster._t_mean_v.size()-braggInfo.first-1;
      braggInfo = std::pair<size_t,double>(newIdx,braggInfo.second);
    }
    // index of hit at which the bragg peak is found
    size_t braggIdx = braggInfo.first;
    // amount of charge (in hit area) at the bragg peak)
    double braggQ   = braggInfo.second;
    // distance between the bragg peak and the end of the track
    double braggDistToEnd = abs( cluster._s_v[braggIdx] - cluster._s_v[cluster._s_v.size()-1] );

    // get the point that is 15 cm away from the Bragg Peak
    size_t MIPendIdx = _calc.GetMIPendPos(cluster._s_v,braggInfo.first,15);
    // if the track isn't long enough, return false (we won't consider
    // this track)
    if (MIPendIdx == michel::kINVALID_SIZE)
      return false;

    // make a vector of the MIP region only
    std::vector<double> muonMIPdS = std::vector<double>(cluster._s_v.begin(),cluster._s_v.begin()+MIPendIdx);
    std::vector<double> muonMIPdQ = std::vector<double>(cluster._t_mean_v.begin(),cluster._t_mean_v.begin()+MIPendIdx);

    // get the median and RMS of the MIP ADC values
    // make a copy 'cause we don't want the original sorted
    auto muonMIPdQsorted = muonMIPdQ;
    double MIPmedian = _calc.GetMedian(muonMIPdQsorted);
    double MIPrms    = _calc.stdev(muonMIPdQ);

    // get the hit indices for MIP hits within 1 RMS of median
    std::vector<size_t> MIPindices = _calc.GetMIPindices(muonMIPdQ,MIPmedian,MIPrms,1.);
    // get the vectors for these indices
    std::vector<double> MIPdS = _calc.GetSubVector(muonMIPdS,MIPindices);
    std::vector<double> MIPdQ = _calc.GetSubVector(muonMIPdQ,MIPindices);
    // get the linear fit values for these vectors
    auto fit = _calc.GetLinearFit(MIPdS,MIPdQ);
    double MIPs = fit.first;
    double MIPm = fit.second;

    // corrected RMS (taking into account only truly MIP hits)
    double MIPrms_corr = _calc.GetRms(MIPdS,MIPdQ,MIPs,MIPm);

    // calcualte bragg amplitude/area
    double braggExpected = MIPm + MIPs * braggIdx;
    double braggAmp = braggQ - braggExpected;
    double braggArea = _calc.GetBraggArea(cluster._s_v,cluster._t_mean_v,MIPendIdx,braggIdx,MIPm,MIPs);

    // now decide what to do...
    // if area < minimum value -> no michel here
    if (braggArea < _minBraggArea)
      return false;

    cluster._boundary = braggBoundary;

    return true;
  }



}
#endif
