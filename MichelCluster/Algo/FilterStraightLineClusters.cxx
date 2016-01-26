#ifndef MICHELCLUSTER_FILTERSTRAIGHTLINECLUSTERS_CXX
#define MICHELCLUSTER_FILTERSTRAIGHTLINECLUSTERS_CXX

#include "FilterStraightLineClusters.h"
#include <sstream>
namespace michel{

  FilterStraightLineClusters::FilterStraightLineClusters(){

    _name = "FilterStraightLineClusters";
    _min_rms = 0.;

  }
  
  void FilterStraightLineClusters::EventReset()
  {}

  bool FilterStraightLineClusters::FilterCluster(MichelCluster& cluster){

    // calculate the RMS for all hits in time / wire separately
    double avg_time = 0;
    double rms_time = 0;
    double avg_wire = 0;
    double rms_wire = 0;
    
    // calculate average
    for (auto const& h : cluster._all_hits){
      avg_time += h._t;
      avg_wire += h._t;
    }
    avg_time /= cluster._all_hits.size();
    avg_wire /= cluster._all_hits.size();

    // calculate rms
    for (auto const& h : cluster._all_hits){
      rms_time += (h._t - avg_time) * (h._t - avg_time);
      rms_wire += (h._w - avg_wire) * (h._w - avg_wire);
    }
    
    rms_time = rms_time/(cluster._all_hits.size()+1);
    rms_wire = rms_wire/(cluster._all_hits.size()+1);

    if ( (rms_time < (_min_rms*_min_rms)) or
	 (rms_wire < (_min_rms*_min_rms)) )
	 return false;
    
    return false;
  }
  
}

#endif
