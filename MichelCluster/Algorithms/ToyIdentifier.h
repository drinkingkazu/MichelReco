/**
 * \file ToyIdentifier.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ToyIdentifier
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_TOYIDENTIFIER_H
#define MICHELCLUSTER_TOYIDENTIFIER_H

#include "Fmwk/BaseAlgIdentifier.h"
namespace michel {
  /**
     \class ToyIdentifier
  */
  class ToyIdentifier : public BaseAlgIdentifier {
    
  public:
    
    /// Default constructor
    ToyIdentifier(){}
    
    /// Default destructor
    ~ToyIdentifier(){}

    /// Event resetter
    void EventReset();

    /// Toy Michel identifier
    Michel Identify(const MichelCluster& cluster);
    
  };
}

#endif
/** @} */ // end of doxygen group 

