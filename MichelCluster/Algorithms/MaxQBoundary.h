/**
 * \file MaxQBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class MaxQBoundary
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_MAXQBOUNDARY_H
#define MICHELCLUSTER_MAXQBOUNDARY_H

#include "Fmwk/BaseAlgBoundary.h"
namespace michel {
  /**
     \class MaxQBoundary
  */
  class MaxQBoundary : public BaseAlgBoundary {
    
  public:
    
    /// Default constructor
    MaxQBoundary(){}
    
    /// Default destructor
    ~MaxQBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(const MichelCluster& cluster);
    
  private:
    size_t find_max(const std::vector<double>& data);

  };
}

#endif
/** @} */ // end of doxygen group 

