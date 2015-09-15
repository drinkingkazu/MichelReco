#ifndef MICHELCLUSTER_REQUIRECOVARIANCEDIP_CXX
#define MICHELCLUSTER_REQUIRECOVARIANCEDIP_CXX

#include "RequireCovarianceDip.h"
#include "Fmwk/MichelException.h"
//#include "ClusterVectorCalculator.h"
#include <cmath>
#include <sstream>

namespace michel {
  
  void RequireCovarianceDip::EventReset()
  {}
  
  bool RequireCovarianceDip::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  { 

    /// May be later I will extend this to N number of dips but it's
    /// probably overkill
    
    if(cluster._chi2_v.size() < cluster._hits.size()) {
      Print(msg::kEXCEPTION,this->Name(),"Covariance vector size less than num hits, run CalcTruncated");
      throw MichelException();
    }

    const auto& covariance = cluster._chi2_v;
    
    /// Loop over covariance array in both directions, if you see something...
    /// say something says MTA on the subway, lets do the same here
    /// look for a SINGLE dip in covariance below cutoff
    
    double cutoff           = _covariance_dip_cutoff;
    bool   in_low_reg       = false;
    bool   been_in_low_reg  = false;
    bool   been_in_high_reg = false;
      
    /// loop one direction only for now, I don't
    /// yet see the benefit to looping the other direction...
    
    for(unsigned int i = 0; i < covariance.size(); ++i) {

      double r = fabs(covariance[i]);

      if(been_in_low_reg && in_low_reg)
	return false;
      
      //starting beginning of low region
      if(r < cutoff) {
	in_low_reg = true;
	continue;
      }

      // am I still in low_region
      if(r < cutoff && in_low_reg) {
	continue;
      }

      //maybe i'm not in low region
      if(r >= cutoff) {
	//but what if I was (in_low_region true upon previous loop), then we been
	if(in_low_reg) {
	  been_in_low_reg = true;
	}
	in_low_reg = false;
	been_in_high_reg = true;
	continue;
      }
      
    }
    
    if(in_low_reg && been_in_high_reg) been_in_low_reg = true;
    
    if(!been_in_low_reg)
      return false;
    
    return true; // Yes there was a single dip in linearity
  }
  
}
#endif
