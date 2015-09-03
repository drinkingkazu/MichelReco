#ifndef MICHEL_CUTONFIDUCIALVOLUME_CXX
#define MICHEL_CUTONFIDUCIALVOLUME_CXX

#include "CutOnFiducialVolume.h"
#include <sstream>

namespace michel {

CutOnFiducialVolume::CutOnFiducialVolume()
{
  _name    = "CutOnFiducialVolume";
  _excluded_wire_ranges.clear();
  _excluded_time_ranges.clear();

  //Default buffer size set to 3cm.
  _wire_buffer_size = 3.;
  _time_buffer_size = 3.;

}

bool CutOnFiducialVolume::ProcessCluster(MichelCluster& cluster,
    const std::vector<HitPt>& hits)
{

  //If user never set ranges to exclude, of course do not exclude any clusters.
  if (!_excluded_wire_ranges.size() && !_excluded_time_ranges.size()) return true;

  //If for some reason the michel doesn't have any hits, of course exclude it
  //(this shouldn't ever happen)
  if (!cluster._michel.size()) return false;

  //Loop over the michel hits and make sure none are within the excluded regions
  for (auto const& ihit : cluster._michel) {
    //For this hit, make sure it is not in any of the excluded wire regions
    for (auto const& _ex_w : _excluded_wire_ranges )
      if (ihit._w > (_ex_w.first - _wire_buffer_size) &&
          ihit._w < (_ex_w.second + _wire_buffer_size) )
        return false;
    //For this hit, make sure it is not in any of the excluded time regions
    for (auto const& _ex_t : _excluded_time_ranges )
      if (ihit._w > (_ex_t.first - _time_buffer_size) &&
          ihit._w < (_ex_t.second + _time_buffer_size) )
        return false;
      
  }

  return true;
}

}

#endif
