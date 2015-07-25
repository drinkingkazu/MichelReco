#ifndef MICHELCLUSTER_MICHEL_CXX
#define MICHELCLUSTER_MICHEL_CXX

#include "Michel.h"
#include <iostream>
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


  void Michel::Dump() {
    std::cout << "\n\t -- Start Michel --\n";
    std::cout << "\tenergy : " << _charge << " length: " << _length << "\n";
    std::cout << "\tat position : (" << _start._w << "," << _start._t << ")" << std::endl;
    std::cout << "\twith number of hits " << this->size() << std::endl;
    std::cout << "\t -- End Michel --\n";
  }
}

#endif
