/**
 * \file BaseMichelAlgo.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseMichelAlgo
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEMICHELALGO_H
#define BASEMICHELALGO_H

#include "MichelTypes.h"
namespace michel {
  /**
     \class BaseMichelAlgo
  */
  class BaseMichelAlgo{
    
  public:
    
    /// Default constructor
    BaseMichelAlgo()
    { _verbosity = msg::kNORMAL; }
    
    /// Default destructor
    virtual ~BaseMichelAlgo(){}

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level)
    { _verbosity = level; }
    
    /// Event-wise reset function
    virtual void EventReset() = 0;

  protected:

    /// Verbosity level
    msg::MSGLevel_t _verbosity;
    
  };
}

#endif
/** @} */ // end of doxygen group 

