/**
 * \file BaseAlgIdentifier.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgIdentifier
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGIDENTIFIER_H
#define BASEALGIDENTIFIER_H

#include "BaseMichelAlgo.h"
#include "MichelCluster.h"
namespace michel {
  /**
     \class BaseAlgIdentifier
     User defined class BaseAlgIdentifier ... these comments are used to generate
     doxygen documentation!
  */
  class BaseAlgIdentifier : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    BaseAlgIdentifier(){}
    
    /// Default destructor
    virtual ~BaseAlgIdentifier(){}

    virtual Michel Identify(const MichelCluster& cluster) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

