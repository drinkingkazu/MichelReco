#ifndef MICHELCLUSTERTYPES_H
#define MICHELCLUSTERTYPES_H

#include <cstdlib>

namespace michel {

  /// Type def for index address
  typedef size_t HitIdx_t;

  /// Type def for hit ID
  typedef size_t HitID_t;

  /// Type def for cluster ID
  typedef size_t ClusterID_t;

  namespace msg {
    /// Verbosity message level
    enum MSGLevel_t {
      kDEBUG,
      kINFO,
      kNORMAL,
      kWARNING,
      kERROR
    };
  }

  ///
  enum AlgoType_t {
    kClusterMerger,
    kBoundaryFinder,
    kMichelID,
    kMichelCluster,
    kMIDFilter,
    kAlgoTypeMax
  };

  struct EventID {
    int run;
    int subrun;
    int event;
  };
  
}
#endif
