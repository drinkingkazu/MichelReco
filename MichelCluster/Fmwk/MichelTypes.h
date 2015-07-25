#ifndef MICHELCLUSTERTYPES_H
#define MICHELCLUSTERTYPES_H

#include <cstdlib>

namespace michel {

  /// Type def for index address
  typedef size_t HitIdx_t;

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
    kMichelFinder,
    kAlgoTypeMax
  };
  
}
#endif
