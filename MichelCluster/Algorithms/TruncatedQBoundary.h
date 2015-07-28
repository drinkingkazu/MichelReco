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

namespace michel {
  /**
     \class TruncatedQBoundary
  */
  class TruncatedQBoundary : public BaseAlgBoundary {
    
  public:
    
    /// Default constructor
    TruncatedQBoundary(){}
    
    /// Default destructor
    ~TruncatedQBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    HitIdx_t Boundary(MichelCluster& cluster);

    
  private:
    unsigned int nCk( unsigned int n, unsigned int k );

    double coeff(double k, double N);
    
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

