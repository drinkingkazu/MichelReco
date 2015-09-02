/**
 * \file TSpectrumBoundary.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class TSpectrumBoundary
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELCLUSTER_TSPECTRUMBOUNDARY_H
#define MICHELCLUSTER_TSPECTRUMBOUNDARY_H

#include <algorithm>

#include "Fmwk/BaseMichelAlgo.h"

#include "TSpectrum.h"
#include "TH1D.h"
#include "TPolyMarker.h"
#include "TList.h"

namespace michel {
  /**
     \class TSpectrumBoundary
  */
  class TSpectrumBoundary : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    TSpectrumBoundary(){_name="TSpectrumBoundary";}
    
    /// Default destructor
    ~TSpectrumBoundary(){}

    /// Event resetter
    void EventReset();
    
    /// A function to identify a michel's boundary point
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

    
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


    size_t get_tspectrum_max(const std::vector<double>& X, 
			     const std::vector<double>& Y);

    size_t find_max(const std::vector<double>& data);
    size_t find_min(const std::vector<double>& data);
    
  };
}

#endif
/** @} */ // end of doxygen group 

