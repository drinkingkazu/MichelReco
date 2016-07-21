/**
 * \file ClusterPhotons.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ClusterPhotons
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CLUSTERPHOTONS_H
#define CLUSTERPHOTONS_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"
#include <map>
#include <sstream>

namespace michel {
  /**
     \class ClusterPhotons
     User defined class ClusterPhotons ... these comments are used to generate
     doxygen documentation!
  */
  class ClusterPhotons : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    ClusterPhotons();
    
    /// Default destructor
    ~ClusterPhotons(){};

    /// Event re-setter
    void EventReset(){};

    /**
     * @brief Cut on michel clusters that have high average charge per hit
     */
    bool ProcessCluster(MichelCluster& michel,
			const std::vector<HitPt>& hits);

    /**
     *@brief set the maximum distance between hits for them to be in the same cluster
     */
    void SetMaxDist(double d) { _d_max = d; }

  private:

    // max distance between hits for them to be in same cluster
    double _d_max;

    // map hit index -> cluster index
    std::map<size_t,size_t> _hit_cluster_map;
    
  };
}

#endif
/** @} */ // end of doxygen group 

