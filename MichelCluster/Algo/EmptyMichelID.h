/**
 * \file EmptyMichelID.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class EmptyMichelID
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_EMPTYMICHELID_H
#define MICHELCLUSTER_EMPTYMICHELID_H

#include "Fmwk/BaseAlgMichelID.h"
namespace michel {
/**
   \class EmptyMichelID
*/
class EmptyMichelID : public BaseAlgMichelID {

public:

  /// Default constructor
  EmptyMichelID() {}

  /// Default destructor
  ~EmptyMichelID() {}

  /// Event resetter
  void EventReset() {};

  /// Toy Michel identifier
  Michel Identify(const MichelCluster& cluster, bool& forward);

};
}

#endif
/** @} */ // end of doxygen group

