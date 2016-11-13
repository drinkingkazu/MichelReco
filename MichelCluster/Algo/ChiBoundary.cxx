#ifndef MICHELCLUSTER_CHIBOUNDARY_CXX
#define MICHELCLUSTER_CHIBOUNDARY_CXX

#include "ChiBoundary.h"
#include "TGraphErrors.h"
#include "TF1.h"


namespace michel {

  void ChiBoundary::EventReset()
  {}
    
  bool ChiBoundary::ProcessCluster(MichelCluster& cluster,
				   const std::vector<HitPt>& hits) {

    if (hits.size() == 0) return false;

    //stolen from run_michel
    float  _chi2_rise = 5;
    float _chi2_fall = 6;
    float  _chi2_thresh = 0;

    //make vector of chi in cluster
    auto chi = do_chi(cluster, 15); //hardcoded window size
    
    //find max chi
    auto the_chi_max_peaks   = find_max_pos(chi,
					    true,
					    15,
					    1.0,
					    _chi2_rise,
					    _chi2_fall,
					    _chi2_thresh);

    HitIdx_t max_chi = kINVALID_SIZE;

    //return index of max, or invalid if none found
    if (the_chi_max_peaks.size() > 0 ) {
      max_chi = the_chi_max_peaks.at(0); //replace with max
    }

    cluster._boundary = cluster._ordered_pts[max_chi];
    return true;
  }
  
  
  std::vector<double> ChiBoundary::do_chi(const MichelCluster& cluster, size_t window_size){
    
    std::vector<double> chi;

    auto chi_data = get_windows(cluster._ordered_pts,window_size);
    
    chi.reserve(chi_data.size());
    std::vector<double> x   (window_size,0.);
    std::vector<double> y   (window_size,0.);
    std::vector<double> xerr(window_size,0.3);
    std::vector<double> yerr(window_size,1.0);
    
    TGraphErrors graph(window_size,&x[0],&y[0],&xerr[0],&yerr[0]);

    TF1 tf("aho","[0] + [1]*x");
  
    for(const auto& cd : chi_data) {
      if(cd.size() != window_size){
	chi.push_back(0.0);
	continue;
      }
    
      for(size_t i = 0; i < window_size; ++i) {
	graph.SetPoint( i,
			cluster._hits[cluster._ordered_pts[i]]._w,
			cluster._hits[cluster._ordered_pts[i]]._t);
      }
      
      tf.SetParameter(0,1);
      tf.SetParameter(1,1);
      
      graph.Fit(&tf,"F 0 N Q");
      
      double amin = tf.GetChisquare()/(window_size - 1);
      chi.push_back(amin);
      
    }
    
    return chi;


  }
  
  //get_windows
  
  
  template<typename T>
  std::vector<std::vector<T> > ChiBoundary::get_windows(const std::vector<T>& the_thing,
							const size_t window_size)
  {
  
    std::vector<std::vector<T> > data;
			  
    auto w = window_size + 2;
    w = (unsigned int)((w - 1)/2);
    auto num = the_thing.size();

    data.reserve(num);
  
    for(size_t i = 1; i <= num; ++i) {
      std::vector<T> inner;
      inner.reserve(20);
      if(i < w) {
	for(size_t j = 0; j < 2 * (i%w) - 1; ++j)
	  inner.push_back(the_thing[j]);
      }else if (i > num - w + 1){
	for(size_t j = num - 2*((num - i)%w)-1 ; j < num; ++j)
	  inner.push_back(the_thing[j]);
      }else{
	for(size_t j = i - w; j < i + w - 1; ++j)
	  inner.push_back(the_thing[j]);
      }
      data.emplace_back(inner);
    }

    return data;
  
  }

  //find_max_pos
  std::vector<HitIdx_t> ChiBoundary::find_max_pos(const std::vector<double> chi,
						  bool forward,
						  size_t window,
						  float cutoff,
						  float rise_edge, 
						  float fall_edge, 
						  float threshold){
    std::vector<HitIdx_t> result;
    auto ped_info = PedEstimate(chi,forward, window, cutoff);
    float ped_mean = ped_info.first;
    float ped_rms = ped_info.second;
    bool found_pulse = false; 
    HitIdx_t t = 0;
  
    HitIdx_t size = chi.size();
  
    while (  t < size ) {
      if(chi[t] > (ped_mean  + rise_edge * ped_rms)  && 
	 chi[t] > (threshold + ped_mean) &&
	 !found_pulse)
	found_pulse = true;
      
      if(found_pulse) {
	HitIdx_t  t_end = t;

	while(1) {
	  if(t_end == size - 1) {
	    ++t_end;
	    goto END;
	  }
	  
	  if(chi[t_end]  <= (ped_mean + fall_edge * ped_rms) &&
	     (chi[t_end] <= threshold + ped_mean))
	    break;
	  else 
	    ++t_end;
	}
	result.push_back(find_max_peak(chi, t, t_end));	
	
      END:
	while(t < t_end) ++t; //secretly increases t...
	
      }
      ++t;
      found_pulse = false;
    }
    
    return result;
  }

  //find max peak
  HitIdx_t ChiBoundary::find_max_peak(const std::vector<double>& chi,
				      HitIdx_t istart,
				      HitIdx_t iend)
  {
    auto the_max = double{0.0};
    HitIdx_t cl = 4096;
  
    for(HitIdx_t i = istart; i < iend; ++i) {
      if(chi[i] > the_max) { the_max = chi[i]; cl = i; }
    }
  
    return cl;
  }


  //ped estimator
  std::pair<float,float> ChiBoundary::PedEstimate(const std::vector<double>& data, bool start, HitIdx_t window, HitIdx_t cutoff) {
    float mean = 0;
    float rms = 0;
    HitIdx_t n = data.size();
    HitIdx_t k = 0;
    //number of points to consider in calculation;
    //int window = 10;
  
    bool below = false;
    //need minimum number to calculate
    if (n >= window) {
      if (start == true){
	while (below == false && k+window < n){
	  auto mean_rms = getrms(data, k, k+window, window);
	  mean = mean_rms.first;
	  rms  = mean_rms.second;
	
	  if (rms < cutoff){
	    below = true;
	  }
	  k++;
	}
      }
      else{
	k = n-1;
	while (below == false && k - window > 0){
	  auto mean_rms =  getrms(data, k-window, k, window);
	  mean = mean_rms.first;
	  rms  = mean_rms.second;
	  if (rms < cutoff){
	    below = true;
	  }
	  k--;
	}
      }
    }
  
    //returns <mean, rms>, or 0,0 if nothing in vector or bad index or not below cutoff
    return std::pair<float,float>(mean,rms);
  }

  std::pair<float,float> ChiBoundary::getrms (const std::vector<double>& data, HitIdx_t k, HitIdx_t m, HitIdx_t window) 
  {
    float mean = 0;
    float rms = 0;
    for (size_t i  = k; i < m; i++) mean += data.at(i);
  
    mean = mean/window;

    for (size_t i= k; i < m; i++){
      float diff = data.at(i) - mean;
      rms += diff*diff;
    }
  
    rms = sqrt(rms/window);
    return  std::pair<float,float>(mean,rms);
  }

  
}
#endif
