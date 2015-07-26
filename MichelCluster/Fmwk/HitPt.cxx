#ifndef MICHELCLUSTER_HITPT_CXX
#define MICHELCLUSTER_HITPT_CXX

#include "HitPt.h"
#include <cmath>
namespace michel {

  HitPt::HitPt(const double q,
	       const double w,
	       const double t,
	       const size_t id)
    : _id ( id )
    , _q  ( q  )
    , _w  ( w  )
    , _t  ( t  )
  {}
  
  double HitPt::SqDist(const HitPt& h) const
  { return (pow(_w - h._w,2) + pow(_t - h._t,2)); }
  
  double HitPt::Dist(const HitPt& h) const
  { return sqrt(SqDist(h)); }
}
#endif
