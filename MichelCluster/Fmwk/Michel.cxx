#ifndef MICHELCLUSTER_MICHEL_CXX
#define MICHELCLUSTER_MICHEL_CXX

#include "Michel.h"
#include "UtilFunc.h"
#include <iostream>
#include <sstream>
namespace michel {

  Michel::Michel()
    : _charge ( kINVALID_DOUBLE )
    , _energy ( kINVALID_DOUBLE )
    , _length ( kINVALID_DOUBLE )
    , _start  ()
  {}
  
  Michel::Michel(const double charge,
		 const double energy,
		 const double length, 
		 const HitPt& start)
    : _charge ( charge )
    , _energy ( energy )
    , _length ( length )
    , _start  ( start  )
  {}

  std::string Michel::Diff(const Michel& rhs) const
  {
    std::string msg;
    if(!IsSame((*this),rhs)) 
      msg += "    Michel hit list changed\n";

    if(_start != rhs._start)
      msg += "    HitPt _start changed\n";

    if(_charge != rhs._charge)
      msg += "    double _charge changed\n";

    if(_length != rhs._length)
      msg += "    double _length changed\n";

    if(_energy != rhs._energy)
      msg += "    double _energy changed\n";

    if(_charge != rhs._charge)
      msg += "    double _charge changed\n";
    return msg;
  }

  void Michel::Dump() const {
    std::stringstream ss;
    ss << "\n\t -- Start Michel --" << std::endl
       << "\tenergy : " << _energy << " length: " << _length << " charge " << _charge << "" << std::endl
       << "\tat position : (" << _start._w << "," << _start._t << ")" << std::endl
       << "\twith number of hits " << this->size() << std::endl
       << "\t -- End Michel --";
    Print(msg::kNORMAL,__FUNCTION__,ss.str());
  }
}

#endif
