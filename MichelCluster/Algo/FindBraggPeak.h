/**
 * \file FindBraggPeak.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class FindBraggPeak
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef FINDBRAGGPEAK_H
#define FINDBRAGGPEAK_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
  /**
     \class FindBraggPeak
  */
  class FindBraggPeak : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    FindBraggPeak();
    
    /// Default destructor
    ~FindBraggPeak(){}

    /// Event re-setter
    void EventReset();

    /// Re-cluster michel electrons w/ un-used hits
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// setter for minimum Bragg Area allowed
    void SetMinBraggArea(double a) { _minBraggArea = a; }

  private:

    /// minimum bragg area allowed for this cluster
    /// to be considered a stopping muon
    double _minBraggArea;
    
  };
}

#endif
/** @} */ // end of doxygen group 

