/**
 * \file EdgeMerger.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class EdgeMerger
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_EDGEMERGER_H
#define MICHELCLUSTER_EDGEMERGER_H

#include "Fmwk/BaseAlgBinaryMerger.h"
namespace michel {
  /**
     \class EdgeMerger
  */
  class EdgeMerger : public BaseAlgBinaryMerger {
  
  public:
    
    /// Default constructor
    EdgeMerger() : _edge_dist(6.0) {}
    
    /// Default destructor
    ~EdgeMerger(){}

    /// Event reset
    void EventReset();

    /// Merge function to assign a pair-wise score for a decision making
    bool Merge(const MichelCluster& a, const MichelCluster& b);

    /// Priority function assign a priority ordering for a merging function to be called
    double Priority(const MichelCluster& cluster);
    
    /// setter function _edge_distance
    void SetEdgeDistance(double d) { _edge_dist = d; }
    
  private:

    /// Check if two clusters are touching 
    bool Touching (const MichelCluster& lhs,
		   const MichelCluster& rhs,
		   const double min_dist) const;
    
    /// Edge distance between two clusters to define "touching"
    double _edge_dist;
  };
  
}

#endif
/** @} */ // end of doxygen group 

