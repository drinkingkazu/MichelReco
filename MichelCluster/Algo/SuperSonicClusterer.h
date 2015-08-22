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

#include "Fmwk/BaseAlgMichelCluster.h"

namespace michel {
  /**
     \class SuperSonicClusterer
  */
  class SuperSonicClusterer : public BaseAlgMichelCluster {
    
  public:
    
    /// Default constructor
    SuperSonicClusterer(){}
    
    /// Default destructor
    ~SuperSonicClusterer(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    void Cluster(Michel& cluster,
		 const std::vector<HitPt>& hits);
    
  };
}

#endif
/** @} */ // end of doxygen group 

