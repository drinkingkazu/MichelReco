/**
 * \file CutOnMuonLinearity.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CutOnMuonLinearity
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CUTONMUONLINEARITY_H
#define CUTONMUONLINEARITY_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class CutOnMuonLinearity
     User defined class CutOnMuonLinearity ... these comments are used to generate
     doxygen documentation!
  */
  class CutOnMuonLinearity : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    CutOnMuonLinearity();
    
    /// Default destructor
    ~CutOnMuonLinearity(){};

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
     *@brief set minimum chi allowed for a hit to be considered in the straight part of muon
     */
    void SetChiMin(double c) { _chi_min = c; }

    /**
     *@brief set minimum fraction of hits in cluster that can be used to get slope (if less then don't use algo)
     */
    void SetFracMinHits(double f) { _frac_min_hits = f; }

  private:

    // a whole bunch of settable parameters

    // minimum chi value allowed for a hit to be considered
    // in the straight section of a muon
    double _chi_min;

    // minimum fraction of all hits in the cluster
    // that need to have a chi above a certain value
    // i.e. a measure of the muon's overall linearity
    double _frac_min_hits;

  };
}

#endif
/** @} */ // end of doxygen group 

