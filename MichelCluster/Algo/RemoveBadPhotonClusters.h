/**
 * \file RemoveBadPhotonClusters.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RemoveBadPhotonClusters
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef REMOVEBADPHOTONCLUSTERS_H
#define REMOVEBADPHOTONCLUSTERS_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class RemoveBadPhotonClusters
     User defined class RemoveBadPhotonClusters ... these comments are used to generate
     doxygen documentation!
  */
  class RemoveBadPhotonClusters : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RemoveBadPhotonClusters();
    
    /// Default destructor
    ~RemoveBadPhotonClusters(){};

    /// Event re-setter
    void EventReset(){};

    /**
     * @brief Use MichelCluster and surrounding hits to decide if this is really a michel
     * @input MichelCluster michel : the currently reconstructed michel object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    bool ProcessCluster(MichelCluster& michel,
			const std::vector<HitPt>& hits);

    /**
     *@brief set the maximum allowed photon cluster linearity
     */
    void SetMaxLinearity(double l) { _max_linearity = l; }

  private:

    // maximum linearity for photon to be considered ok
    double _max_linearity;
    
  };
}

#endif
/** @} */ // end of doxygen group 

