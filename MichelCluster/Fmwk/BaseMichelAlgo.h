/**
 * \file BaseMichelAlgo.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseMichelAlgo
 *
 * @author kazuhiro + david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEMICHELALGO_H
#define BASEMICHELALGO_H

#include "MichelTypes.h"
#include "MichelCluster.h"

namespace michel {
  /**
     \class BaseMichelAlgo
  */
  class BaseMichelAlgo{
    
  public:
    
    /// Default constructor
    BaseMichelAlgo()
      { _verbosity = msg::kNORMAL; _name = "BaseMichelAlgo"; }
    
    /// Default destructor
    virtual ~BaseMichelAlgo(){}

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level)
    { _verbosity = level; }
    
    /// Event-wise reset function
    virtual void EventReset() = 0;

    /**
       @brief Algorithm function to edit MichelCluster object
       @details This function is the main tool for all algorithms
       in this framework. It takes an input (editable) MichelCluster
       and a const reference to a hit-list containing all hits in the
       event. The algorithm should do something to improve the
       MichelCluster or decide that this in fact is not a Michel.
     */
    virtual bool ProcessCluster(MichelCluster& cluster,
				const std::vector<HitPt>& hits) = 0;

    /**
     * @brief Return name of algorithm
     */
    const std::string Name() const { return _name.c_str(); }

  protected:

    /// Verbosity level
    msg::MSGLevel_t _verbosity;

    /// Name for algorithm
    std::string _name;
    
  };
}

#endif
/** @} */ // end of doxygen group 

