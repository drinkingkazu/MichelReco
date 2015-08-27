/**
 * \file EmptyMichelCluster.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class EmptyMichelCluster
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef EMPTYMICHELCLUSTER_H
#define EMPTYMICHELCLUSTER_H

#include "Fmwk/BaseAlgMichelCluster.h"

namespace michel {
/**
   \class EmptyMichelCluster
*/
class EmptyMichelCluster : public BaseAlgMichelCluster {

public:

  /// Default constructor
  EmptyMichelCluster() {}

  /// Default destructor
  ~EmptyMichelCluster() {}

  /// Event re-setter
  void EventReset() {};

  /// Re-cluster michel electrons w/ un-used hits
  void Cluster(Michel& cluster,
               const std::vector<HitPt>& hits);

};
}

#endif
/** @} */ // end of doxygen group

