#ifndef REMOVEFAKEPMTSIGNALS_CXX
#define REMOVEFAKEPMTSIGNALS_CXX

#include "RemoveFakePMTSignals.h"
#include "Fmwk/ClusterVectorCalculator.h"

#include <sstream>

namespace michel{

  RemoveFakePMTSignals::RemoveFakePMTSignals()
  {
    _name    = "RemoveFakePMTSignals";
    _maxErrorTime = 0.;
  }

  bool RemoveFakePMTSignals::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)

  {

    if (hits.size() == 0) return false;

    ClusterVectorCalculator _clusterCalc;

    // check the time of all hits in the michel
    // if almost all of them have the exact same time-tick
    // they are probably due to these weird PMT signals
    // that cross onto the wires
    // remove the michel
    auto const& michel = cluster._michel;

    // get a vector of the times in the hits
    std::vector<double> michel_times;
    for (auto const& h : michel)
      michel_times.push_back(h._t);

    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "michel Times: ";
      for (auto const& t : michel_times) { ss << t << ", "; }
      Print(msg::kINFO,this->Name(),ss.str());
    }

    // get the median time
    double medianT = _clusterCalc.GetMedian(michel_times);
    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "median is " << medianT;
      Print(msg::kINFO,this->Name(),ss.str());
    }

    // what fraction of hits have the median time value?
    int count = 0;
    for (auto const& t : michel_times)
      if (medianT == t) { count += 1; }

    double frac_median = (double)count / (double)michel_times.size();
    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << frac_median << " is frac. of hits that have the median time.";
      Print(msg::kINFO,this->Name(),ss.str());
    }

    // calculate the RMS time-value
    // if too low (like really low)
    // then probably a perfectly horizontal PMT signal
    double avgT = 0;
    double rmsT = 0;
    for (auto const& t : michel_times)
      avgT += t;
    avgT /= (double)michel_times.size();

    for (auto const& t : michel_times)
      rmsT += (t-avgT)*(t-avgT);
    rmsT = sqrt(rmsT/(double)michel_times.size());

    double uncertainty = rmsT/sqrt((double)michel_times.size());

    if (_verbosity <= msg::kINFO){
      std::stringstream ss;
      ss << "RMS on time is " << rmsT << "\tUncertainty on Mean is: " << uncertainty;
      Print(msg::kINFO,this->Name(),ss.str());
    }

    if (uncertainty < _maxErrorTime){
      if (_verbosity <= msg::kINFO){
	std::stringstream ss;
	ss << "Removing this event!";
	Print(msg::kINFO,this->Name(),ss.str());
      }
      return false;
    }

    // if rms is low enough!
    return true;
  }

}

#endif
