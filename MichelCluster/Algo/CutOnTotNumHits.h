/**
 * \file CutOnTotNumHits.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CutOnTotNumHits
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CUTONTOTNUMHITS_H
#define CUTONTOTNUMHITS_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class CutOnTotNumHits
     User defined class CutOnTotNumHits ... these comments are used to generate
     doxygen documentation!
  */
  class CutOnTotNumHits : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    CutOnTotNumHits();
    
    /// Default destructor
    ~CutOnTotNumHits(){};

    /// Event re-setter
    void EventReset(){};

    /**
     * @brief Cut on the total number of hits in the event
     * @detail If the number of hits is less than some threshold
     then this cluster should be ignored. This algorithm should be run
     immediately after the ClusterMerging stage, to get rid of any
     surviving small clusters such as ugly clumps of hits in the event
     or delta-rays that should not be scanned.
     * @input MichelCluster michel : the currently reconstructed michel object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    bool ProcessCluster(MichelCluster& michel,
			const std::vector<HitPt>& hits);

    /**
     *@brief set the minimum number of hits the entire cluster should have
     */
    void SetMinTotHits(int n) { _min_num_hits = n; }

  private:

    // minimum number of hits a cluster can have
    double _min_num_hits;

  };
}

#endif
/** @} */ // end of doxygen group 

