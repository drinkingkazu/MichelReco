/**
 * \file ToyBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ToyBoundary
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_TOYBOUNDARY_H
#define MICHELCLUSTER_TOYBOUNDARY_H

#include "Fmwk/BaseAlgBoundary.h"
namespace michel {
  /**
     \class ToyBoundary
  */
  class ToyBoundary : public BaseAlgBoundary {
    
  public:
    
    /// Default constructor
    ToyBoundary(){}
    
    /// Default destructor
    ~ToyBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(MichelCluster& cluster);
    
  };
}

#endif
/** @} */ // end of doxygen group 

