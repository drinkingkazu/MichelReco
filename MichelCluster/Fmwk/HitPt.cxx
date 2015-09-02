#ifndef MICHELCLUSTER_HITPT_CXX
#define MICHELCLUSTER_HITPT_CXX

#include "HitPt.h"
#include <cmath>
#include <sstream>

namespace michel {

  HitPt::HitPt(const double  q,
	       const double  w,
	       const double  t,
	       const HitID_t id,
	       const int     pl)
    : _id ( id )
    , _pl ( pl )
    , _q  ( q  )
    , _w  ( w  )
    , _t  ( t  )
  {}
  
  double HitPt::SqDist(const HitPt& h) const
  { return (pow(_w - h._w,2) + pow(_t - h._t,2)); }
  
  double HitPt::Dist(const HitPt& h) const
  { return sqrt(SqDist(h)); }

  std::string HitPt::Print() const
  {
    std::stringstream ss;
    ss << " ID : " << _id << " ... "
       << " PL = " << _pl << " : "
       << " Q = " << _q << " : "
       << " W = " << _w << " : "
       << " T = " << _t << std::endl;
    return ss.str();
  }
}
#endif
