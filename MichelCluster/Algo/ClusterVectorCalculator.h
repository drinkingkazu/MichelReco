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
#include "Fmwk/ColorPrint.h"
namespace michel {
  /**
     \class ClusterVectorCalculator
  */
  class ClusterVectorCalculator : public ColorPrint{
    
  public:
    
    /// Default constructor
    ClusterVectorCalculator(){}
    
    /// Default destructor
    ~ClusterVectorCalculator(){}

    std::vector<double> calc_smooth_derive(const std::vector<double>& _dist,
					   const std::vector<double>& tmeans, 

					   const int s) const;

    std::vector<double> calc_smooth_mean  (const MichelCluster& cluster,
					   const double _n_window_size,
					   const int window_cutoff,
					   const double p_above) const;
    
    double do_smooth_derive(const std::vector<double>& f,
			    const std::vector<double>& x,
			    int N) const;

    /// calculate slope from least squares at point...
    std::vector<double> calc_slope(const std::vector<::michel::HitPt>& hits,
				   const int _n_window_size) const;


    unsigned int nCk( unsigned int n, unsigned int k ) const;

    double coeff(double k, double N) const;


    template<typename S>
    S calc_mean(const std::vector<S>& data) const;
    
    template<typename W>
    void cut(std::vector<W>& data, double frac) const;
    
    template<typename T>
    std::vector<std::vector<T> > get_windows(const std::vector<T>& the_thing,
					     const int window_size) const;



    std::vector<double> calc_covariance(const std::vector<::michel::HitPt>& hits, 
					const int _n_window_size) const;

    
    double cov (const std::vector<double>& data1,
		const std::vector<double>& data2) const;

    double stdev(const std::vector<double>& data) const;
    double mean (const std::vector<double>& data) const;

    size_t find_max(const std::vector<double>& data) const;
    size_t find_min(const std::vector<double>& data) const;

    // function to calculate linear fit to list of points
    std::pair<double,double> GetLinearFit(const std::vector<michel::HitPt>& pts) const;

    // function to calculate linear fit to list of points
    std::pair<double,double> GetLinearFit(const std::vector<double>& x,
					  const std::vector<double>& y);

    double GetPerpendicularDistance(const michel::HitPt& h,
				    const double& slope,
				    const double& intercept) const;

    double GetMedian(std::vector<double>& v) const;


    // get the position marked as the end of the MIP region
    // this is 15 cm before the bragg peak is reached
    size_t GetMIPendPos(const std::vector<double>& v,
			const size_t& max,
			const double distAsked) const;

    // get indices for hits considered MIP (within some alpha*RMS of the median charge)
    std::vector<size_t> GetMIPindices(const std::vector<double>& dQ,
				      const double& median,
				      const double& rms,
				      const double& alpha) const;

    // get a subvector containing the entries for the list of indices provided
    std::vector<double> GetSubVector(const std::vector<double> v,
				     const std::vector<size_t> idx) const;

    // given (x,y) points calculate the RMS w.r.t. a linear
    // equation y = s * x + b
    // basically get rms of residuals
    double GetRms(const std::vector<double>& MIPdS,
		  const std::vector<double>& MIPdQ,
		  const double& slope,
		  const double& intercept) const;

    // get max index and value
    std::pair<size_t,double> GetMaxIndex(const std::vector<double>& v) const;

    // function to calculate Bragg Peak area under the curve (taking into account best fit to MIP slope)
    double GetBraggArea(const std::vector<double>& dS, const std::vector<double>& dQ,
			const size_t& MIPendIdx, const size_t& braggIdx,
			const double& MIPm, const double& MIPs) const;
    
  };
}

#endif
/** @} */ // end of doxygen group 

