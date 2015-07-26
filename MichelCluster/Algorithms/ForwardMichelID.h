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
    ForwardMichelID(){}
    
    /// Default destructor
    ~ForwardMichelID(){}

    /// Event resetter
    void EventReset();

    /// Toy Michel identifier
    Michel Identify(const MichelCluster& cluster);

  private:
    double determine_length(const MichelCluster& c,
			    const std::vector<size_t> ordered_pts_idx);
    
    bool determine_forward(const MichelCluster& cluster,
			   const double n_cutoff,
			   const double c_cutoff,
			   const double w_cutoff,
			   bool&  forward);
    
  };
}

#endif
/** @} */ // end of doxygen group 

