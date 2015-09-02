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

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class StepSuperSonicCluster
  */
  class StepSuperSonicCluster : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    StepSuperSonicCluster() : _max_step(3.6) {}
    
    /// Default destructor
    ~StepSuperSonicCluster(){_name="StepSuperSonicCluster";}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
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

