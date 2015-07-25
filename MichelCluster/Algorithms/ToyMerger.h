/**
 * \file ToyMerger.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ToyMerger
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_TOYMERGER_H
#define MICHELCLUSTER_TOYMERGER_H

#include "Fmwk/BaseAlgMerger.h"
namespace michel {
  /**
     \class ToyMerger
  */
  class ToyMerger : public BaseAlgMerger {
  
  public:
    
    /// Default constructor
    ToyMerger(){}
    
    /// Default destructor
    ~ToyMerger(){}

    /// Event reset
    void EventReset();

    /// Merge clusters
    MichelClusterArray Merge(const MichelClusterArray& input_v);
  };
    
}

#endif
/** @} */ // end of doxygen group 

