/**
 * \file ClusterVectorCalculator.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class ClusterVectorCalculator
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_CLUSTERVECTORCALCULATOR_H
#define MICHELCLUSTER_CLUSTERVECTORCALCULATOR_H

#include <algorithm>
#include "Fmwk/MichelCluster.h"

namespace michel {
  /**
     \class ClusterVectorCalculator
  */
  class ClusterVectorCalculator {
    
  public:
    
    /// Default constructor
    ClusterVectorCalculator(){}
    
    /// Default destructor
    ~ClusterVectorCalculator(){}

    std::vector<double> calc_smooth_derive(const std::vector<double>& _dist,
					   const std::vector<double>& tmeans, 

					   const int s);
    std::vector<double> calc_smooth_mean  (const MichelCluster& cluster,
					   const double _n_window_size,
					   const int window_cutoff,
					   const double p_above);
    
    double do_smooth_derive(const std::vector<double>& f,
			    const std::vector<double>& x,
			    int N);

    /// calculate slope from least squares at point...
    std::vector<double> calc_slope(const std::vector<::michel::HitPt>& hits,
				   const int _n_window_size);


    unsigned int nCk( unsigned int n, unsigned int k );

    double coeff(double k, double N);


    template<typename S>
    S calc_mean(const std::vector<S>& data);
    
    template<typename W>
    void cut(std::vector<W>& data, double frac);
    
    template<typename T>
    std::vector<std::vector<T> > get_windows(const std::vector<T>& the_thing,
					     const int window_size);



    std::vector<double> calc_covariance(const std::vector<::michel::HitPt>& hits, 
					const int _n_window_size);
    
    double cov (const std::vector<double>& data1,
		const std::vector<double>& data2);

    double stdev(const std::vector<double>& data);
    double mean (const std::vector<double>& data);

    size_t find_max(const std::vector<double>& data);
    size_t find_min(const std::vector<double>& data);
    
  };
}

#endif
/** @} */ // end of doxygen group 

