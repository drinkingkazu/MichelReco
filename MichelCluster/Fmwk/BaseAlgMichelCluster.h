/**
 * \file BaseAlgMichelCluster.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgMichelCluster
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGMICHELCLUSTER_H
#define BASEALGMICHELCLUSTER_H

#include "BaseMichelAlgo.h"
#include "MichelCluster.h"
namespace michel {
  /**
     \class BaseAlgMichelCluster
     User defined class BaseAlgMichelCluster ... these comments are used to generate
     doxygen documentation!
  */
  class BaseAlgMichelCluster : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    BaseAlgMichelCluster(){}
    
    /// Default destructor
    virtual ~BaseAlgMichelCluster(){}

    virtual void Cluster(Michel& michel,
			 const std::vector<HitPt>& hits) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

