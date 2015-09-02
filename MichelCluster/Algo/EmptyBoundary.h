/**
 * \file EmptyBoundary.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class EmptyBoundary
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_EMPTYBOUNDARY_H
#define MICHELCLUSTER_EMPTYBOUNDARY_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
/**
   \class EmptyBoundary
*/
  class EmptyBoundary : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    EmptyBoundary() {}
    
    /// Default destructor
    ~EmptyBoundary() {}
    
    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    void ProcessCluster( MichelCluster& cluster,
			 const std::vector<HitPt>& hits);
  };
}

#endif
/** @} */ // end of doxygen group

