/**
 * \file RemoveBraggPeakHits.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RemoveBraggPeakHits
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef REMOVEBRAGGPEAKHITS_H
#define REMOVEBRAGGPEAKHITS_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class RemoveBraggPeakHits
  */
  class RemoveBraggPeakHits : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RemoveBraggPeakHits();
    
    /// Default destructor
    ~RemoveBraggPeakHits(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// set max radius within which to search for potential Bragg hits
    void SetMaxRadius(double r) { _sq_radius = r*r; }

    /// set maximum charge allowed ( factor represents x times above mean after which hit is
    /// considered for exclusion )
    void SetChargeFactor(double f) { _f = f; }

  private:

    ///  if hit Q > mean Q * _f -> hit is candidate for removal
    double _f;

    /// maximum distance^2 from michel boundary point
    /// where to search for hits to possibly be removed
    double _sq_radius;

  };
}

#endif
/** @} */ // end of doxygen group 

