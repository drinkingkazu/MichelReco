/**
 * \file BaseAlgBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgBoundary
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGBOUNDARY_H
#define BASEALGBOUNDARY_H

#include "BaseMichelAlgo.h"
#include "MichelCluster.h"
namespace michel {
  /**
     \class BaseAlgBoundary
  */
  class BaseAlgBoundary : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    BaseAlgBoundary(){}
    
    /// Default destructor
    virtual ~BaseAlgBoundary(){}

    /// A function to identify a michel's boundary point
    virtual HitIdx_t Boundary(const MichelCluster& cluster) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

