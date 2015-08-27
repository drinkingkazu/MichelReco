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

#include "Fmwk/BaseAlgBoundary.h"

namespace michel {
/**
   \class EmptyBoundary
*/
class EmptyBoundary : public BaseAlgBoundary {

public:

  /// Default constructor
  EmptyBoundary() {}

  /// Default destructor
  ~EmptyBoundary() {}

  /// Event resetter
  void EventReset();

  /// A function to identify a michel's boundary point
  HitIdx_t Boundary( MichelCluster& cluster);

};
}

#endif
/** @} */ // end of doxygen group

