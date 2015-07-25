/**
 * \file ToyMichelID.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ToyMichelID
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_TOYMICHELID_H
#define MICHELCLUSTER_TOYMICHELID_H

#include "Fmwk/BaseAlgMichelID.h"
namespace michel {
  /**
     \class ToyMichelID
  */
  class ToyMichelID : public BaseAlgMichelID {
    
  public:
    
    /// Default constructor
    ToyMichelID(){}
    
    /// Default destructor
    ~ToyMichelID(){}

    /// Event resetter
    void EventReset();

    /// Toy Michel identifier
    Michel Identify(const MichelCluster& cluster);
    
  };
}

#endif
/** @} */ // end of doxygen group 

