/**
 * \file SuperSonicClusterer.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class SuperSonicClusterer
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef SUPERSONICCLUSTERER_H
#define SUPERSONICCLUSTERER_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class SuperSonicClusterer
  */
  class SuperSonicClusterer : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    SuperSonicClusterer();
    
    /// Default destructor
    ~SuperSonicClusterer(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// Setter for merge-till-converge
    void SetMergeTillConverge(bool on) { _merge_till_converge = on; }

    /// set max radius within which to search for michel hits
    void SetMaxRadius(double r) { _max_radius = r; }

  private:

    /// Merge 'till converge flag
    bool _merge_till_converge;

    /// max radius within to search for potential hits to
    /// add to michel
    double _max_radius;
    
  };
}

#endif
/** @} */ // end of doxygen group 

