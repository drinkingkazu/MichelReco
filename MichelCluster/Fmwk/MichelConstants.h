#ifndef MICHELCLUSTERCONSTANTS_H
#define MICHELCLUSTERCONSTANTS_H

#include <climits>
#include <limits>

namespace michel {

  static const double kINVALID_DOUBLE = std::numeric_limits<double>::max();

  static const double kMAX_DOUBLE = std::numeric_limits<double>::max();
  
  static const double kMIN_DOUBLE = std::numeric_limits<double>::min();
  
  static const std::size_t kINVALID_SIZE = std::numeric_limits<std::size_t>::max();

  static const int kINVALID_INT = std::numeric_limits<int>::max();

}
#endif
