/**
 * \file TruncatedQBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class TruncatedQBoundary
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_TRUNCATEDQBOUNDARY_H
#define MICHELCLUSTER_TRUNCATEDQBOUNDARY_H

#include "Fmwk/BaseAlgBoundary.h"
#include <algorithm>

namespace michel {
  /**
     \class TruncatedQBoundary
  */
  class TruncatedQBoundary : public BaseAlgBoundary {
    
  public:
    
    /// Default constructor
    TruncatedQBoundary(){_maxDistance=20;}
    
    /// Default destructor
    ~TruncatedQBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(MichelCluster& cluster);

    /// setter function for max distance between min in dQds
    /// and max in dQ
    void SetMaxDistanceTruncatedPeaks(int n) { _maxDistance = n; }

    
  private:

    /// retuns binomial coefficient
    unsigned int nCk( unsigned int n, unsigned int k );

    /// utility function for do_smooth_derive function, 
    /// return coefficient for special difference scheme
    double coeff(double k, double N);
    
    /// return vector of smoothed derivative
    std::vector<double> calc_smooth_derive(const std::vector<double>& _dist,
					   const std::vector<double>& tmeans, 
					   const int s);

    /// return vector of truncated mean
    std::vector<double> calc_smooth_mean  (const MichelCluster& cluster,
					   const double _n_window_size,
					   const int window_cutoff,
					   const double p_above);
    
    
    /// calculate smooth differivative at point
    double do_smooth_derive(const std::vector<double>& f,
			    const std::vector<double>& x,
			    int N);
    
    
    /// calculate slope from least squares at point...
    std::vector<double> calc_slope(const std::vector<::michel::HitPt>& hits,
				   const int _n_window_size);
      
    /// redundant function, same as mean function below
    template<typename S>
    S calc_mean(const std::vector<S>& data);
    
    /// delete upper frac*len(data) in data
    template<typename W>
    void cut(std::vector<W>& data, double frac);

    /// return vector of the_thing, binned into length window_size
    template<typename T>
    std::vector<std::vector<T> > get_windows(const std::vector<T>& the_thing,
					     const int window_size);



    /// return vector of covariance, in window _n_window_size
    std::vector<double> calc_covariance(const std::vector<::michel::HitPt>& hits, 
					const int _n_window_size);
    

    /// return covariance between two data sets
    double cov (const std::vector<double>& data1,
		const std::vector<double>& data2);

    /// return stdev
    double stdev(const std::vector<double>& data);

    /// return average
    double mean (const std::vector<double>& data);
    
    /// return index of maximum value in data
    size_t find_max(const std::vector<double>& data);

    /// return index of minimum value in data
    size_t find_min(const std::vector<double>& data);

    /// number of hits that sets max allowed distance between
    /// min in trunc. dQdx and max in trunc. dQ
    /// if the min and max are more than this number of hits
    /// apart -> do not create a michel
    int _maxDistance;
    
  };
}

#endif
/** @} */ // end of doxygen group 

