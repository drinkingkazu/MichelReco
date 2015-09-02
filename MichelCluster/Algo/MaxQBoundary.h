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

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class MaxQBoundary
  */
  class MaxQBoundary : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    MaxQBoundary(){_name="MaxQBoundary";}
    
    /// Default destructor
    ~MaxQBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

  private:
    size_t find_max(const std::vector<double>& data);

  };
}

#endif
/** @} */ // end of doxygen group 

