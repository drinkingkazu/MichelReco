/**
 * \file RequireSlopeSignFlip.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RequireSlopeSignFlip
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_REQUIRESLOPESIGNFLIP_H
#define MICHELCLUSTER_REQUIRESLOPESIGNFLIP_H

#include "Fmwk/BaseMichelAlgo.h"
#include <algorithm>

namespace michel {
  /**
     \class RequireSlopeSignFlip
  */
  class RequireSlopeSignFlip : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RequireSlopeSignFlip() {_name="RequireSlopeSignFlip";}
    
    /// Default destructor
    ~RequireSlopeSignFlip(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);
    
  private:
    
    /// return sign of val (-1/1 and sometimes 0)
    int sign(double val);
    
  };
  
  
}

#endif
/** @} */ // end of doxygen group 

