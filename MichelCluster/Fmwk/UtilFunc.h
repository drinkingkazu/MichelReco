#ifndef MICHELRECO_UTILFUNC_H
#define MICHELRECO_UTILFUNC_H

#include <vector>
#include <cstddef>

namespace michel {

  template <class T>
  bool IsSame(const std::vector<T>& lhs, const std::vector<T> rhs)
  {
    if(lhs.size() != rhs.size()) return false;
    for(size_t i=0; i<lhs.size(); ++i) 

      if(lhs[i] != rhs[i]) return false;

    return true;
  }

}

#endif
