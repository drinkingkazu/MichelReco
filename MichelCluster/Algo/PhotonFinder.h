/**
 * \file PhotonFinder.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class PhotonFinder
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef PHOTONFINDER_H
#define PHOTONFINDER_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class PhotonFinder
  */
  class PhotonFinder : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    PhotonFinder();
    
    /// Default destructor
    ~PhotonFinder(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// set max radius within which to search for michel hits
    void SetMaxRadius(double r) { _max_radius = r; }

    /// set minimum dot-product
    void SetMinDotProduct(double d) { _min_dot = d; }

  private:

    /// minimum dot-product
    double _min_dot;

    /// maximum distance from michel cluster for a hit to even
    /// be considered
    double _max_radius;

  };
}

#endif
/** @} */ // end of doxygen group 

