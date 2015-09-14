/**
 * \file CalcTruncated.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CalcTruncated
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_REQUIRECOVARIANCEDIP_H
#define MICHELCLUSTER_REQUIRECOVARIANCEDIP_H

#include "Fmwk/BaseMichelAlgo.h"
#include <algorithm>

namespace michel {
  /**
     \class RequireCovarianceDip
  */
  class RequireCovarianceDip : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    RequireCovarianceDip()
      :
      _covariance_dip_cutoff(0.9)
    {_name="RequireCovarianceDip";}
      
      /// Default destructor
      ~RequireCovarianceDip(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);
    
    /// setter function for strength SINGLE dip in covariance must be
    /// 0.9 is default
    void SetCovarianceDipCutoff(double d)    { _covariance_dip_cutoff = d; }
    
  private:
    
    /// single dip in covariance strenght, must be below this value to
    /// be considered a dip
    double _covariance_dip_cutoff;
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

