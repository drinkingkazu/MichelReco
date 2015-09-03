#ifndef MICHELCLUSTERTYPES_H
#define MICHELCLUSTERTYPES_H

#include <cstdlib>
#include <string>
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
      kERROR,
      kEXCEPTION,
      kMSG_TYPE_MAX
    };
    const std::string kColorPrefix[kMSG_TYPE_MAX] =
      {
	"\033[94m", ///< blue ... DEBUG
	"\033[92m", ///< green ... INFO
	"\033[95m", ///< magenta ... NORMAL
	"\033[93m", ///< yellow ... WARNING
	"\033[91m", ///< red ... ERROR
	"\033[5;1;33;41m" ///< red with yellow background ... CRITICAL
      };
    ///< Color coding of message

    const std::string kStringPrefix[kMSG_TYPE_MAX] =
      {
	"     [DEBUG]  ", ///< DEBUG message prefix
	"      [INFO]  ", ///< INFO message prefix
	"    [NORMAL]  ", ///< NORMAL message prefix
	"   [WARNING]  ", ///< WARNING message prefix
	"     [ERROR]  ", ///< ERROR message prefix
	" [EXCEPTION]  "  ///< CRITICAL message prefix
      };
    ///< Prefix of message
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
