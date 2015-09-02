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

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class RadiusMichelCluster
  */
  class RadiusMichelCluster : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
  RadiusMichelCluster() : _min_radius(10.0), _max_hit_charge(500.0) {_name="RadiusMichelCluster";}
    
    /// Default destructor
    ~RadiusMichelCluster(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);
    
    /// Setter for minimum radius
    void SetMinRadius(double r) { _min_radius = r; }

    /// Setter for maximum allowed single hit charge
    void SetMaxHitCharge(double q) { _max_hit_charge = q; }
    
  private:
    
    /// Minimum radius from michel start point.
    double _min_radius;

    /// Maximum allowable charge for a single hit to be added to a Michel.
    /// This is set to try to avoid including the last (high-charge) hit from the muon bragg peak.
    double _max_hit_charge;

  };
}

#endif
/** @} */ // end of doxygen group 

