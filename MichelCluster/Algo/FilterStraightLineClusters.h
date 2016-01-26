/**
 * \file FilterStraightLineClusters.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class FilterStraightLineClusters
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster
    
    @{*/
#ifndef MICHELCLUSTER_FILTERSTRAIGHTLINECLUSTERS_H
#define MICHELCLUSTER_FILTERSTRAIGHTLINECLUSTERS_H

#include "Fmwk/BaseAlgFilter.h"
namespace michel {

  /**
     \class FilterStraightLineClusters
  */
  class FilterStraightLineClusters : public BaseAlgFilter {
  
  public:
    
    /// Default constructor
    FilterStraightLineClusters();
    
    /// Default destructor
    ~FilterStraightLineClusters(){}

    /// Event reset
    void EventReset();

    /// Merge function to assign a pair-wise score for a decision making
    bool FilterCluster(const MichelCluster& cluster);

    /// set minimum RMS allowed for a cluster (to avoid straight lines)
    void setMinRMS(double m) { _min_rms = m; }

  private:

    double _min_rms;

  };
  
}

#endif
/** @} */ // end of doxygen group 

