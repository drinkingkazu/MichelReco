#ifndef RADIUSMICHELCLUSTER_CXX
#define RADIUSMICHELCLUSTER_CXX

#include "RadiusMichelCluster.h"

namespace michel {

void RadiusMichelCluster::EventReset()
{}

void RadiusMichelCluster::ProcessCluster(MichelCluster& cluster,
					 const std::vector<HitPt>& hits) {

  auto& michel = cluster._michel;

  /// This michel was bogus when it came in, don't cluster further
  if (michel.size() == 0) return;

  double radius = michel._length;

  /// Set radius to min rad of michel isn't long enough
  if (radius < _min_radius) radius = _min_radius;

  if (_verbosity <=  msg::kINFO) {
    std::cout << "\n\n\tRadius clustering"
              << "\n\twe use a radius of: " << radius
              << "\n\tto cluster the remaining "
              << "\n\thits.size() : " << hits.size()
              << "\n\tinput\n\n";
  }


  auto michel_start = michel._start;

  for (const auto& thishit : hits)
    if (michel_start.SqDist(thishit) <= radius * radius)
      if (thishit._q < _max_hit_charge)
        michel.push_back(thishit);

  michel._charge = 0;

  for (const auto& michel_hit : michel)
    michel._charge += michel_hit._q;

}



}
#endif
