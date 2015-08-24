/**
 * \file ForwardMichelID.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ForwardMichelID
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_FORWARDMICHELID_H
#define MICHELCLUSTER_FORWARDMICHELID_H

#include "Fmwk/BaseAlgMichelID.h"
namespace michel {
  /**
     \class ForwardMichelID
  */
  class ForwardMichelID : public BaseAlgMichelID {
    
  public:
    
    /// Default constructor
    ForwardMichelID(){_maxHits=0;}
    
    /// Default destructor
    ~ForwardMichelID(){}

    /// Event resetter
    void EventReset();

    /// Toy Michel identifier
    Michel Identify(const MichelCluster& cluster, bool& forward);

    /// Setter for max number of hits allowed
    void SetMaxMichelHits(size_t n) { _maxHits = n; std::cout << "max hits: " << _maxHits << std::endl; }

  private:
    double determine_length(const MichelCluster& c,
			    const std::vector<size_t>& ordered_pts_idx);
    
    bool determine_forward(const MichelCluster& cluster,
			   const double n_cutoff,
			   const double c_cutoff,
			   const double w_cutoff,
			   bool&  forward);

    /// variable to determine max number of hits allowed in michel cluster
    size_t _maxHits;
    
  };
}

#endif
/** @} */ // end of doxygen group 

