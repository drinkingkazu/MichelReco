/**
 * \file BaseAlgMerger.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgMerger
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGMERGER_H
#define BASEALGMERGER_H

#include "BaseMichelAlgo.h"
#include "MichelCluster.h"
namespace michel {
  /**
     \class BaseAlgMerger
  */
  class BaseAlgMerger : public BaseMichelAlgo {
  
  public:
    
    /// Default constructor
    BaseAlgMerger() : BaseMichelAlgo() {}
    
    /// Default destructor
    virtual ~BaseAlgMerger(){}

    /// Function to be implemented by children classes
    virtual MichelClusterArray Merge(const MichelClusterArray& input_v) = 0;

  };
    
}

#endif
/** @} */ // end of doxygen group 

