/**
 * \file StepAroundCluster.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class StepAroundCluster
 *
 * @author KaZuHirO
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef STEPAROUNDCLUSTER_H
#define STEPAROUNDCLUSTER_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class StepAroundCluster
  */
  class StepAroundCluster : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
  StepAroundCluster() : _max_step(3.6) {_name="StepAroundCluster.h";}
    
    /// Default destructor
    ~StepAroundCluster(){}

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

