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

#include "MichelCluster.h"

namespace michel {
  /**
     \class BaseAlgMerger
  */
  class BaseAlgMerger{
  
  public:
    
    /// Default constructor
    BaseAlgMerger()
      { _verbosity = msg::kNORMAL; _name = "BaseAlgMerger"; }
    
    /// Default destructor
    virtual ~BaseAlgMerger(){}

    /// Event-wise reset function
    virtual void EventReset() = 0;

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level)
    { _verbosity = level; }

    /// Function to be implemented by children classes
    virtual MichelClusterArray Merge(const MichelClusterArray& input_v) = 0;

    /**
     * @brief Return name of algorithm
     */
    std::string Name() { return _name; }

  protected:

    /// Verbosity level
    msg::MSGLevel_t _verbosity;

    /// Name for algorithm
    std::string _name;

  };
    
}

#endif
/** @} */ // end of doxygen group 

