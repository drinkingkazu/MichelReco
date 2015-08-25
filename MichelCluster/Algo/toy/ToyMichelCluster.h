/**
 * \file ToyMichelCluster.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ToyMichelCluster
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef TOYMICHELCLUSTER_H
#define TOYMICHELCLUSTER_H

#include "Fmwk/BaseAlgMichelCluster.h"
namespace michel {
  /**
     \class ToyMichelCluster
  */
  class ToyMichelCluster : public BaseAlgMichelCluster {
    
  public:
    
    /// Default constructor
    ToyMichelCluster(){}
    
    /// Default destructor
    ~ToyMichelCluster(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    void Cluster(Michel& michel,
		 const std::vector<HitPt>& hits);
    
  };
}

#endif
/** @} */ // end of doxygen group 

