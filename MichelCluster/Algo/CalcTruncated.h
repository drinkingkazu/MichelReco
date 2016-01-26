/**
 * \file CalcTruncated.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CalcTruncated
 *
 * @author vic
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_CALCTRUNCATED_H
#define MICHELCLUSTER_CALCTRUNCATED_H

#include "Fmwk/BaseMichelAlgo.h"
#include <algorithm>
#include <TStopwatch.h>

namespace michel {
  /**
     \class CalcTruncated
  */
  class CalcTruncated : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    CalcTruncated();
    
    /// Default destructor
    ~CalcTruncated(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    /// setter function for window size for calculating covariance and slope
    /// this in principle (and currently is) different from truncated Q window size
    void SetCovarianceWindowSize(int s)      { _covariance_window = s; }
    
    /// setter function for window size used to calculated truncated Q
    void SetTruncatedQWindowSize(int s)      { _n_window_size     = s; }
    
    /// setter function for fraction of points to keep in truncated Q calc
    void SetPAbove(double p)                 { _p_above = p; }
    
    // setter function for truncated Q window cutoff, windows less than
    // this size have no truncation
    void SetMinWindowSize(int w)             { _window_cutoff = w; }
    
    /// setter function for mitigating edge effect, this number of points on the
    /// edge set equal to truncatedQ[_edgefix], maybe 0 is best I don't know
    void SetEdgeEffectFix(int e)             { _edgefix = e; }

    /// report resource usage
    void Report();
    
  private:

    /// Sliding window size for calculating covariance/slope window
    int _covariance_window;
    
    /// Sliding window size for calculating truncated Q window
    int _n_window_size;
    
    /// Fraction of points to keep in truncated Q calculation, for example
    /// a value of 0.25 means we keep 25% of the lowest points in the specified
    /// window size, then do mean
    double _p_above;
    
    /// Minimum window size to apply truncation, if window size is smaller than this
    /// then we keep all the points 
    int  _window_cutoff;
    
    /// Try to correct edge effect by setting this number of points equal
    /// on both edges
    int _edgefix;
    
    // stop-watch for time-profiling
    TStopwatch _watch;
    double _t0, _t1, _t2, _t3;
    int    _n0, _n1, _n2, _n3;
    
  };

  
}

#endif
/** @} */ // end of doxygen group 

