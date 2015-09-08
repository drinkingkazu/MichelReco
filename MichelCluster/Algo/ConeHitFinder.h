/**
 * \file ConeHitFinder.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ConeHitFinder
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CONEHITFINDER_H
#define CONEHITFINDER_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class ConeHitFinder
  */
  class ConeHitFinder : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    ConeHitFinder();
    
    /// Default destructor
    ~ConeHitFinder(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// set max radius within which to search for michel hits
    void SetMaxRadius(double r) { _max_radius = r; }

    /// set maximum perpendicular distance to line
    void SetMaxPerpendicularDistance(double d) { _max_perp_distance = d; }

  private:

    /// maximum perpendicular distance to the line
    double _max_perp_distance;

    /// maximum distance from michel cluster for a hit to even
    /// be considered
    double _max_radius;

  };
}

#endif
/** @} */ // end of doxygen group 

