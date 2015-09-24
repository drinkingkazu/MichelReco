/**
 * \file RequireLargeAngle.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RequireLargeAngle
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef REQUIRELARGEANGLE_H
#define REQUIRELARGEANGLE_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class RequireLargeAngle
     User defined class RequireLargeAngle ... these comments are used to generate
     doxygen documentation!
  */
  class RequireLargeAngle : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RequireLargeAngle();
    
    /// Default destructor
    ~RequireLargeAngle(){};

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
     * @brief set the minimum angle between muon and michel clusters
     */
    void SetMinAngle(double a) { _min_angle = a; }

    /**
     * @brief minimum number of michel straight hits
     */
    void SetMinStraightMichelHits(int n) { _min_straight_michel_hits = n; }
    
    /**
     * @brief muon distance to use for slope calculation
     */
    void SetMuonLengthUsed(double d) { _muon_length_used = d; }
    
  private:

    // minimum angle between muon and michel clusters
    double _min_angle;

    // minimum number of straight michel hits
    int _min_straight_michel_hits;

    // distance of muon to use for slope calculation
    // starting at the michel boundary
    // if the muon track is smaller than this distance
    // use the full track
    double _muon_length_used;

  };
}

#endif
/** @} */ // end of doxygen group 

