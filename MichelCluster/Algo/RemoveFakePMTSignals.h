/**
 * \file RemoveFakePMTSignals.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RemoveFakePMTSignals
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef REMOVEFAKEPMTSIGNALS_H
#define REMOVEFAKEPMTSIGNALS_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class RemoveFakePMTSignals
     User defined class RemoveFakePMTSignals ... these comments are used to generate
     doxygen documentation!
  */
  class RemoveFakePMTSignals : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RemoveFakePMTSignals();
    
    /// Default destructor
    ~RemoveFakePMTSignals(){};

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
     *@brief set the maximum allowed RMS on time
     */
    void SetMaxErrorTime(double rms) { _maxErrorTime = rms; }

  private:

    // maximum rms allowed on time axis for this in michel
    double _maxErrorTime;
    
  };
}

#endif
/** @} */ // end of doxygen group 

