#ifndef MICHELCLUSTER_REQUIRESLOPESIGNFLIP_CXX
#define MICHELCLUSTER_REQUIRESLOPESIGNFLIP_CXX

#include "RequireSlopeSignFlip.h"
#include "Fmwk/MichelException.h"
#include <cmath>
#include <sstream>

namespace michel {
  
  void RequireSlopeSignFlip::EventReset()
  {}
  
  bool RequireSlopeSignFlip::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  { 
    if(cluster._dirs_v.size() < cluster._hits.size()) {
      Print(msg::kEXCEPTION,this->Name(),"slope vector size less than num hits, run CalcTruncated");
      throw MichelException();
    }
     
    const auto& slope = cluster._dirs_v;
    
    /// Lets just see if the slope changes sign, change in sign means
    /// something changed directions
    
    int  curr_sign = -1;
    int  prev_sign = -1;
    bool changed_sign = false;
    
    for(unsigned int i = 0; i < slope.size(); ++i)
      {
	curr_sign = sign(slope[i]);
	
	if(i == 0) { prev_sign = curr_sign; continue; }
	
	if(curr_sign == prev_sign) {
	  continue;
	}
	else{
	  changed_sign = true;
	  break;
	}
      }
  
    if(!changed_sign) return false;    

    return true;
  }
  
  
  int RequireSlopeSignFlip::sign(double val)
  {
    if (val > 0) return  1;
    if (val < 0) return -1;
    return 0;
  }
  

}
#endif
