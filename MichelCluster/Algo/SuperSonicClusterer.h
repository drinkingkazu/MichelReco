/**
 * \file SuperSonicClusterer.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class SuperSonicClusterer
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef SUPERSONICCLUSTERER_H
#define SUPERSONICCLUSTERER_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class SuperSonicClusterer
  */
  class SuperSonicClusterer : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    SuperSonicClusterer(){_name="SuperSonicClusterer";}
    
    /// Default destructor
    ~SuperSonicClusterer(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);
    
  };
}

#endif
/** @} */ // end of doxygen group 

