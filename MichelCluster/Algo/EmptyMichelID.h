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

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
/**
   \class EmptyMichelID
*/
class EmptyMichelID : public BaseMichelAlgo {

public:

  /// Default constructor
  EmptyMichelID() {_name="EmptyMichelID";}

  /// Default destructor
  ~EmptyMichelID() {}

  /// Event resetter
  void EventReset() {};

  /// Toy Michel identifier
  void ProcessCluster(MichelCluster& cluster,
		      const std::vector<HitPt>& hits);

};
}

#endif
/** @} */ // end of doxygen group

