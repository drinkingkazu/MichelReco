/**
 * \file RadiusMichelCluster.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RadiusMichelCluster
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef RADIUSMICHELCLUSTER_H
#define RADIUSMICHELCLUSTER_H

#include "Fmwk/BaseAlgMichelCluster.h"

namespace michel {
  /**
     \class RadiusMichelCluster
  */
  class RadiusMichelCluster : public BaseAlgMichelCluster {
    
  public:
    
    /// Default constructor
    RadiusMichelCluster() : _min_radius(10.0) {}
    
    /// Default destructor
    ~RadiusMichelCluster(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    void Cluster(Michel& cluster,
		 const std::vector<HitPt>& hits);
    
    /// Setter for minimum radius
    void SetMinRadius(double r) { _min_radius = r; }

    
  private:
    
    /// Minimum radius from michel start point.
    double _min_radius;
    
  };
}

#endif
/** @} */ // end of doxygen group 

