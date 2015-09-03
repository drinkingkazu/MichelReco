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
    StepSuperSonicCluster();
    
    /// Default destructor
    ~StepSuperSonicCluster(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);
      
    /// Setter for minimum radius
    void SetMaxStep(double s) { _max_step = s; }

    /// Setter for merge-till-converge
    void SetMergeTillConverge(bool on) { _merge_till_converge = on; }

    /// set max radius within which to search for michel hits
    void SetMaxRadius(double r) { _max_radius = r; }

    /// hit radius setter
    void SetHitRadius(double r) { _hit_radius = r; }

    /// setter for wether to use hit radius
    void SetUseHitRadius(bool on) { _use_hit_radius = on; }

  private:
    
    /// Minimum radius from michel start point.
    double _max_step;
    
    /// Merge 'till converge flag
    bool _merge_till_converge;

    /// max radius within to search for potential hits to
    /// add to michel
    double _max_radius;

    // hit radius : distance from each michel hit within which
    // to search for new michel hits
    double _hit_radius;

    /// boolean on whether to use fixed hit radius or distance to start
    bool _use_hit_radius;
    
  };
}

#endif
/** @} */ // end of doxygen group 

