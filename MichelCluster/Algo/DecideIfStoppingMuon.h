/**
 * \file DecideIfStoppingMuon.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class DecideIfStoppingMuon
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef DECIDEIFSTOPPINGMUON_H
#define DECIDEIFSTOPPINGMUON_H

#include "Fmwk/BaseAlgMIDFilter.h"
#include "math.h"

namespace michel {
  /**
     \class DecideIfStoppingMuon
     User defined class DecideIfStoppingMuon ... these comments are used to generate
     doxygen documentation!
  */
  class DecideIfStoppingMuon : public BaseAlgMIDFilter {
    
  public:
    
    /// Default constructor
    DecideIfStoppingMuon();
    
    /// Default destructor
    ~DecideIfStoppingMuon(){};

    /// Event re-setter
    void EventReset(){};

    /**
     * @brief Use MichelCluster and surrounding hits to decide if this is really a michel
     * @input MichelCluster michel : the currently reconstructed michel object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    bool IsMichel(const MichelCluster& michel,
		  const std::vector<HitPt>& hits);

    /**
     *@brief set minimum chi allowed for a hit to be considered in the straight part of muon
     */
    void SetChiMin(double c) { _chi_min = c; }

    /**
     *@brief set minimum fraction of hits in cluster that can be used to get slope (if less then don't use algo)
     */
    void SetFracMinHits(double f) { _frac_min_hits = f; }

    /**
     *@brief set max radius within which to check if hits are alligned with muon or not (radius center is michel start)
     */
    void SetHitRadius(double r) { _hit_radius = r; }

    /**
     *@brief set max perpendicular distance from muon-line for a hit to be considered as well aligned with it
     */
    void SetMaxDist(double d) { _max_dist = d; }

    /**
     *@brief minimum number of "bad" hits needed for the michel to be a MID
     */
    void SetMinBadHits(double n) { _min_bad_hits = n; }

  private:

    // a whole bunch of settable parameters

    // minimum chi value allowed for a hit to be considered
    // in the straight section of a muon
    double _chi_min;

    // minimum fraction of all hits in the cluster
    // that needs to be used in order to have a good slope
    double _frac_min_hits;

    // radius centered at the michel start point where to search
    // for hits that can possibly be tagged as aligned with the muon
    double _hit_radius;

    // maximum perpendicular distance to muon-line for a hit
    // to be considered as being compatible with that line
    double _max_dist;

    // minimum number of bad hits needed for the michel to be a MID
    int _min_bad_hits;
    
  };
}

#endif
/** @} */ // end of doxygen group 

