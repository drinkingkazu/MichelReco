/**
 * \file HitPt.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class HitPt
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_HITPT_H
#define MICHELCLUSTER_HITPT_H

#include "MichelConstants.h"
#include "MichelTypes.h"
namespace michel {
  /**
     \class HitPt
     Represents 2D hit point
  */
  class HitPt{
    
  public:
    
    /// Default constructor
    HitPt(const double q = kINVALID_DOUBLE,
	  const double w = kINVALID_DOUBLE,
	  const double t = kINVALID_DOUBLE);
    
    /// Default destructor
    ~HitPt(){}

    double SqDist(const HitPt& h) const;

    double Dist(const HitPt& h) const;

    HitID_t _id; ///< Unique ID
    double  _q;  ///< Charge in ADC*Ticks scale
    double  _w;  ///< Wire in [cm] scale
    double  _t;  ///< Time in [cm] scale
  };
}
#endif
/** @} */ // end of doxygen group 

