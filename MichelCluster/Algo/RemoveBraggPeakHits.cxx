#ifndef REMOVEBRAGGPEAKHITS_CXX
#define REMOVEBRAGGPEAKHITS_CXX

#include "RemoveBraggPeakHits.h"
#include <sstream>
#include <cmath>
#include "Fmwk/ClusterVectorCalculator.h"

namespace michel {

  RemoveBraggPeakHits::RemoveBraggPeakHits()
  {
    _name = "RemoveBraggPeakHits";
    _f      = 2.;
    _sq_radius = 1.; // cm^2
  }

  void RemoveBraggPeakHits::EventReset()
  {}
  
  bool RemoveBraggPeakHits::ProcessCluster(MichelCluster& cluster,
					   const std::vector<HitPt>& hits){

    auto& michel = cluster._michel;

    // How does this algorithm work:
    // given the tagged Bragg peak position in the cluster
    // search for all michel hits that are within a radius _radius
    // from this point
    // any hit within this distance that have a high charge
    // are removed from the Michel cluster
    // as these are likely hits from the Muon Bragg peak
    // that have been mis-clustered
    // "too much charge" = hitQ > meanQ * _f
    // with meanQ the average charge / hit in the Michel cluster hits.
    
    //This michel was bogus when it came in, don't cluster further
    if (michel.size() == 0) return false;


    // get the Michel-Muon boundary position
    double wB = cluster._hits[cluster._boundary]._w;
    double tB = cluster._hits[cluster._boundary]._t;

    // get mean charge of Michel Cluster / hit
    double Qavg = 0.;
    for (auto const& h : michel)
      Qavg += h._q;
    Qavg /= michel.size();

    // loop over Michel hits.
    // if too close to Michel-Muon boundary and too much charge -> remove
    for (size_t n=0; n < michel.size(); n++){
      double w = michel[n]._w;
      double t = michel[n]._t;
      double q = michel[n]._q;
      // if the hit has more than _f*Qavg charge -> consider for removal
      if (q > Qavg * _f)
	continue;
      double d = (w-wB)*(w-wB) + (t-tB)*(t-tB);
      if (d < _sq_radius)
	michel.erase(michel.begin()+n);
    }/// for all hits in the Michel cluster

    return true;
  }


}
#endif
