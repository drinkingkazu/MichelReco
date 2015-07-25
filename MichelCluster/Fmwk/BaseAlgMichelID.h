/**
 * \file BaseAlgMichelID.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgMichelID
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGMICHELID_H
#define BASEALGMICHELID_H

#include "BaseMichelAlgo.h"
#include "MichelCluster.h"
namespace michel {
  /**
     \class BaseAlgMichelID
     User defined class BaseAlgMichelID ... these comments are used to generate
     doxygen documentation!
  */
  class BaseAlgMichelID : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    BaseAlgMichelID(){}
    
    /// Default destructor
    virtual ~BaseAlgMichelID(){}

    virtual Michel Identify(const MichelCluster& cluster) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

