/**
 * \file BaseAlgMIDFilter.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgMIDFilter
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGMIDFILTER_H
#define BASEALGMIDFILTER_H

#include "BaseMichelAlgo.h"
#include "MichelCluster.h"

namespace michel {
  /**
     \class BaseAlgMIDFilter
     User defined class BaseAlgMIDFilter ... these comments are used to generate
     doxygen documentation!
  */
  class BaseAlgMIDFilter : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    BaseAlgMIDFilter(){}
    
    /// Default destructor
    virtual ~BaseAlgMIDFilter(){}

   /**
     * @brief Use MichelCluster and surrounding hits to decide if this is really a michel
     * @input MichelCluster michel : the currently reconstructed michel object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    virtual bool IsMichel(const MichelCluster& michel,
			  const std::vector<HitPt>& hits) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

