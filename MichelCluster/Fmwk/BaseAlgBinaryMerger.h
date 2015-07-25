/**
 * \file BaseAlgBinaryMerger.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgBinaryMerger
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGBINARYMERGER_H
#define BASEALGBINARYMERGER_H

#include "BaseAlgMerger.h"
namespace michel {
  /**
     \class BaseAlgBinaryMerger
  */
  class BaseAlgBinaryMerger : public BaseAlgMerger {
  
  public:
    
    /// Default constructor
    BaseAlgBinaryMerger() : BaseAlgMerger() {}
    
    /// Default destructor
    virtual ~BaseAlgBinaryMerger(){}

    /// BinaryMerger implements specific way of merging by calling underlying Merge and Priority functions
    MichelClusterArray Merge(const MichelClusterArray& input_v);

    /// Merge function to assign a pair-wise score for a decision making (TO BE IMPLEMENTED)
    virtual bool Merge(const MichelCluster& a, const MichelCluster& b) = 0;

    /// Priority function assign a priority ordering for a merging function to be called (TO BE IMPLEMENTED)
    virtual double Priority(const MichelCluster& cluster) = 0;

    /// Recursive merge flag
    void Recursive(bool doit)
    { _recursive = doit; }      

  protected:

    /// Recursive merge flag
    bool _recursive;
  };
    
}

#endif
/** @} */ // end of doxygen group 

