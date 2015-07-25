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
namespace michel {
  /**
     \class HitPt
     User defined class HitPt ... these comments are used to generate
     doxygen documentation!
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

    double _q; ///< Charge in ADC*Ticks scale
    double _w; ///< Wire in [cm] scale
    double _t; ///< Time in [cm] scale
  };
}
#endif
/** @} */ // end of doxygen group 

