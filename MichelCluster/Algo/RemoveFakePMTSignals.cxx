#ifndef REMOVEFAKEPMTSIGNALS_CXX
#define REMOVEFAKEPMTSIGNALS_CXX

#include "RemoveFakePMTSignals.h"
#include <sstream>

namespace michel{

  RemoveFakePMTSignals::RemoveFakePMTSignals()
  {
    _name    = "RemoveFakePMTSignals";
    _maxRMStime = 0.;
  }

  bool RemoveFakePMTSignals::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  {

    // check the time of all hits in the michel
    // if almost all of them have the exact same time-tick
    // they are probably due to these weird PMT signals
    // that cross onto the wires
    // remove the michel
    auto const& michel = cluster._michel;

    // calculate the RMS time-value
    // if too low (like really low)
    // then probably a perfectly horizontal PMT signal
    double avgT = 0;
    double rmsT = 0;
    for (auto const& h : michel)
      avgT += h._t;
    avgT /= (double)michel.size();

    for (auto const& h : michel)
      rmsT += (h._t-avgT)*(h._t-avgT);
    rmsT = sqrt(rmsT/(double)michel.size());

    if (_verbosity <= msg::kINFO)
      std::cout << "RMS on time is " << rmsT << std::endl;

    if (rmsT < _maxRMStime)
      return false;

    // if rms is low enough!
    return true;
  }

}

#endif
