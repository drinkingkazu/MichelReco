/**
 * \file AlwaysNotMID.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class AlwaysNotMID
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef ALWAYSNOTMID_H
#define ALWAYSNOTMID_H

#include "Fmwk/BaseAlgMIDFilter.h"
#include "math.h"

namespace michel {
  /**
     \class AlwaysNotMID
     User defined class AlwaysNotMID ... these comments are used to generate
     doxygen documentation!
  */
  class AlwaysNotMID : public BaseAlgMIDFilter {
    
  public:
    
    /// Default constructor
    AlwaysNotMID();
    
    /// Default destructor
    ~AlwaysNotMID(){};

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

  private:
    
  };
}

#endif
/** @} */ // end of doxygen group 

