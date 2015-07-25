/**
 * \file RadiusMichelCluster.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class RadiusMichelCluster
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef RADIUSMICHELCLUSTER_H
#define RADIUSMICHELCLUSTER_H

#include "Fmwk/BaseAlgMichelCluster.h"
namespace michel {
  /**
     \class RadiusMichelCluster
  */
  class RadiusMichelCluster : public BaseAlgMichelCluster {
    
  public:
    
    /// Default constructor
    RadiusMichelCluster(){}
    
    /// Default destructor
    ~RadiusMichelCluster(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    Michel Cluster(const Michel& cluster,
		   const std::vector<HitPt>& hits);
    
  };
}

#endif
/** @} */ // end of doxygen group 

