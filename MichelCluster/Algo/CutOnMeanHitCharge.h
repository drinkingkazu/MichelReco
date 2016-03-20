/**
 * \file CutOnMeanHitCharge.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CutOnMeanHitCharge
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CUTONMEANHITCHARGE_H
#define CUTONMEANHITCHARGE_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class CutOnMeanHitCharge
     User defined class CutOnMeanHitCharge ... these comments are used to generate
     doxygen documentation!
  */
  class CutOnMeanHitCharge : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    CutOnMeanHitCharge();
    
    /// Default destructor
    ~CutOnMeanHitCharge(){};

    /// Event re-setter
    void EventReset(){};

    /**
     * @brief Cut on michel clusters that have high average charge per hit
     */
    bool ProcessCluster(MichelCluster& michel,
			const std::vector<HitPt>& hits);

    /**
     *@brief set the maximum charge / hit allowed for the Michel cluster
     */
    void SetMaxAvgQ(double q) { _max_qavg = q; }

  private:

    // maximum average charge per hit for Michel clusters
    double _max_qavg;
    
  };
}

#endif
/** @} */ // end of doxygen group 

