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
    _minBraggArea = 1000.;
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

    // Ths algorithm calculates the area of the Bragg Peak area for a cluster

    // 1) find muon's stopping location and select subset of hits belonging to the muon only
    std::vector<double> dS_v;
    std::vector<double> dQ_v;

    auto forward  = cluster._forward;
    auto boundary = cluster._boundary;
    if (forward){
      for (size_t i=0; i < boundary; i++){
	dS_v.push_back(      cluster._s_v[i] );
	dQ_v.push_back( cluster._t_mean_v[i] );
      }
    }
    else{
      for (size_t i = boundary; i < cluster._s_v.size(); i++){
	dS_v.push_back(      cluster._s_v[i] );
	dQ_v.push_back( cluster._t_mean_v[i] );
      }
    }

    // next find the point that is 20 cm from the stopping point.
    // we will consider this to be the start of the bragg peak
    size_t endMIP_idx = 0;
    size_t dS_stop    = dS_v[dS_v.size()-1];
    for (size_t i=0; i < dS_v.size(); i++){
      if ( fabs( dS_v[i] - dS_stop ) > 20 )
	endMIP_idx = i;
    }
    
    // select hits belonging to MIP region only
    std::vector<double> dS_MIP;
    std::vector<double> dQ_MIP;

    for (size_t i=0; i < endMIP_idx; i++){
      dS_MIP.push_back( dS_v[i] );
      dQ_MIP.push_back( dQ_v[i] );
    }


    // calculate a linear fit to the MIP region
    auto fit = _calc.GetLinearFit(dS_MIP,dQ_MIP);
    // slope
    double s = fit.first;
    // intercept
    double m = fit.second;
    
    // calculate the RMS from the fit to the linear region
    double rms = 0.;
    for (size_t i=0; i < dS_MIP.size(); i++){
      double d = (dS_MIP[i] * s + m) - dQ_MIP[i];
      rms += d*d;
    }

    // if RMS == 0 -> no linear fit -> don't trust this track
    // -> return false
    if (rms == 0)
      return false;

    rms = sqrt( rms / dS_MIP.size() );

    // now we proceed backwards from the stopping point
    // and integrate the Bragg Peak area up until
    // the significance of the hit's amplitude
    // is > 3 sigma w.r.t. the expectation for MIPs
    double braggArea = 0.;
    for (size_t i= dQ_v.size()-1; i > endMIP_idx; i--){
      auto ds = dS_v[i];
      auto dq = dQ_v[i];
      auto dq_lin = ds * s + m;
      double sigma = (dq-dq_lin) / rms;
      if (sigma < 3.)
	break;
      else
	braggArea += dq - dq_lin;
    }
    
    // cut on the Bragg Peak area
    if (braggArea < _minBraggArea)
      return false;

    return true;
  }



}
#endif
