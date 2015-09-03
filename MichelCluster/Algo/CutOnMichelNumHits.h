/**
 * \file CutOnMichelNumHits.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CutOnMichelNumHits
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CUTONMICHELNUMHITS_H
#define CUTONMICHELNUMHITS_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class CutOnMichelNumHits
     User defined class CutOnMichelNumHits ... these comments are used to generate
     doxygen documentation!
  */
  class CutOnMichelNumHits : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    CutOnMichelNumHits();
    
    /// Default destructor
    ~CutOnMichelNumHits(){};

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
     *@brief set the minimum number of hits a michel cluster can have
     */
    void SetMinMichelHits(int n) { _min_num_hits = n; }

    /**
     *@brief set the maximum number of hits a michel clsuter can have
     */
    void SetMaxMichelHits(int n) { _max_num_hits = n; }

  private:

    // minimum number of hits a michel cluster can have
    double _min_num_hits;

    // maximum number of hits a michel cluster can have
    double _max_num_hits;
    
  };
}

#endif
/** @} */ // end of doxygen group 

