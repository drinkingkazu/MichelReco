/**
 * \file StepSuperSonicCluster.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class StepSuperSonicCluster
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef STEPSUPERSONICCLUSTER_H
#define STEPSUPERSONICCLUSTER_H

#include "Fmwk/BaseAlgMichelCluster.h"

namespace michel {
  /**
     \class StepSuperSonicCluster
  */
  class StepSuperSonicCluster : public BaseAlgMichelCluster {
    
  public:
    
    /// Default constructor
    StepSuperSonicCluster() : _max_step(3.6) {}
    
    /// Default destructor
    ~StepSuperSonicCluster(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    void Cluster(Michel& cluster,
		 const std::vector<HitPt>& hits);
      
    /// Setter for minimum radius
    void SetMaxStep(double s) { _max_step = s; }

  private:
    
    /// Minimum radius from michel start point.
    double _max_step;
    
  };
}

#endif
/** @} */ // end of doxygen group 

