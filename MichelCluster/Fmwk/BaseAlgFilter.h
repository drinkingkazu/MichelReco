/**
 * \file BaseAlgFilter.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class BaseAlgFilter
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef BASEALGFILTER_H
#define BASEALGFILTER_H

#include "MichelCluster.h"
#include "ColorPrint.h"
namespace michel {
  /**
     \class BaseAlgFilter
  */
  class BaseAlgFilter : public ColorPrint{
  
  public:
    
    /// Default constructor
    BaseAlgFilter()
      { _verbosity = msg::kNORMAL; _name = "BaseAlgFilter"; }
    
    /// Default destructor
    virtual ~BaseAlgFilter(){}

    /// Event-wise reset function
    virtual void EventReset() = 0;

    /**
       @brief Algorithm function that decides if to keep or ignore this cluster
       @details This function will decide, given only an individual MichelCluster's
       information, whether this cluster should be kept and used in downstream analysis
       @input MichelCluster cluster -> the cluster to be exhamined
       @return bool -> whether to keep the event (true) or remove (false)
    */
    
    /// Function to be implemented by children classes
    virtual bool FilterCluster(const MichelCluster& cluster) = 0;

    /**
     * @brief Return name of algorithm
     */
    std::string Name() { return _name; }

  protected:

    /// Name for algorithm
    std::string _name;

  };
    
}

#endif
/** @} */ // end of doxygen group 

