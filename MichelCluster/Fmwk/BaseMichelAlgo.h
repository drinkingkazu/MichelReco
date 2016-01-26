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
#include "ColorPrint.h"
namespace michel {
  /**
     \class BaseMichelAlgo
  */
  class BaseMichelAlgo : public ColorPrint{
    
  public:
    
    /// Default constructor
    BaseMichelAlgo()
      { _verbosity = msg::kNORMAL; _name = "BaseMichelAlgo"; }
    
    /// Default destructor
    virtual ~BaseMichelAlgo(){}

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

    virtual void Report(){}

  protected:

    /// Name for algorithm
    std::string _name;
    
  };
}

#endif
/** @} */ // end of doxygen group 

