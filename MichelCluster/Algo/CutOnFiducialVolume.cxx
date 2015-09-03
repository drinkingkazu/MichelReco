#ifndef MICHEL_CUTONFIDUCIALVOLUME_CXX
#define MICHEL_CUTONFIDUCIALVOLUME_CXX

#include "CutOnFiducialVolume.h"
#include "Fmwk/MichelException.h"
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

  _debug = false;

}

bool CutOnFiducialVolume::ProcessCluster(MichelCluster& cluster,
    const std::vector<HitPt>& hits)
{

  if (_debug)
    std::cout << "CutOnFiducialVolume is excluding "
              << _excluded_wire_ranges.size()
              << " wire ranges, and "
              << _excluded_wire_ranges.size()
              << " time ranges."
              << std::endl;

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


  /// Setter to exclude many wire ranges at once
  void CutOnFiducialVolume::SetExcludedWireRanges(std::vector<double> myranges_min,
						  std::vector<double> myranges_max)
  {
    if (myranges_min.size() != myranges_max.size()){
      std::cout << "\033[93m[ERROR]\033[00m list of min/max wire ranges should be equal in length but are not...you screwed up..."  << std::endl;
      throw MichelException();
  }
    _excluded_wire_ranges.clear();
    for (size_t i=0; i < myranges_min.size(); i++){
      std::pair<double,double> my_pair(myranges_min[i],myranges_max[i]);
      AddExcludedWireRange(my_pair);
    }
  }
  
  /// Setter to exclude many time ranges at once
  void CutOnFiducialVolume::SetExcludedTimeRanges(std::vector<double> myranges_min,
						  std::vector<double> myranges_max)
  {
    if (myranges_min.size() != myranges_max.size()){
      std::cout << "\033[93m[ERROR]\033[00m list of min/max wire ranges should be equal in length but are not...you screwed up..."  << std::endl;
      throw MichelException();
  }
    _excluded_time_ranges.clear();
    for (size_t i=0; i < myranges_min.size(); i++){
      std::pair<double,double> my_pair(myranges_min[i],myranges_max[i]);
      AddExcludedTimeRange(my_pair);
    }
  }
  
}

#endif
