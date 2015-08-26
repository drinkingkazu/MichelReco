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

#include "Fmwk/BaseAlgMichelCluster.h"

namespace michel {
  /**
     \class StepAroundCluster
  */
  class StepAroundCluster : public BaseAlgMichelCluster {
    
  public:
    
    /// Default constructor
    StepAroundCluster() : _max_step(3.6) {}
    
    /// Default destructor
    ~StepAroundCluster(){}

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
