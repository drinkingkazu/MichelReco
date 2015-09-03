#ifndef MICHELCLUSTER_TSPECTRUMBOUNDARY_CXX
#define MICHELCLUSTER_TSPECTRUMBOUNDARY_CXX

#include "TSpectrumBoundary.h"
#include <cmath>
#include <sstream>
namespace michel {
  
  void TSpectrumBoundary::EventReset()
  {}
  
  bool TSpectrumBoundary::ProcessCluster(MichelCluster& cluster,
					 const std::vector<HitPt>& hits)
  { 
    
    std::vector<double> truncated_mean;
    std::vector<double> truncated_dqds;
    
    truncated_mean.reserve(cluster._ordered_pts.size());
    truncated_dqds.reserve(cluster._ordered_pts.size());
    
    //hardcoded for now will become configurable
    double _n_window_size = 15;
    double _p_above       = 0.25;
    int    _window_cutoff = 3;

    //do truncated mean
    truncated_mean = calc_smooth_mean(cluster,
				      _n_window_size,
				      _window_cutoff,
				      _p_above);
    
    for(int i = 0 ; i < 3; ++i) {
      truncated_mean.at(i) = truncated_mean[3];
      truncated_mean.at(truncated_mean.size() - i - 1) = truncated_mean[truncated_mean.size() - 3];
    }
    
    
    
    int s = 3; // must be odd, currently has no setter, sorry that this method has no info on it, ask vic
    truncated_dqds = calc_smooth_derive(cluster._s_v,truncated_mean,s);
    
    //Lets play with truncated mean shaving...
    
    if(_verbosity <= msg::kINFO) {
      std::stringstream ss;
      ss << "\t\tIn TSpectrumBoundary" << std::endl
	 << "\tI have " << truncated_mean.size() << " truncated mean size" << std::endl
	 << "\twith   " << truncated_dqds.size() << " derivative points." << std::endl
	 << "\tMy incoming cluster has " << cluster.GetHits().size() << " hits in it...";
      Print(msg::kINFO,__FUNCTION__,ss.str());
    }
    
    
    
    //auto candidate_loc     = find_max(truncated_mean);
    //auto dqdscandidate_loc = find_min(truncated_dqds); 
    
    //swap empty vector with filled vector
    
    //We need to flip truncated dqds over to find peaks
    std::vector<double> flipped_truncated_dqds;
    flipped_truncated_dqds.reserve(truncated_dqds.size());

    for(auto& dqds : truncated_dqds)
      flipped_truncated_dqds.push_back(-1.0 * dqds);
    
    int candidate_loc     = get_tspectrum_max(cluster._s_v,truncated_mean);
    int dqdscandidate_loc = get_tspectrum_max(cluster._s_v,flipped_truncated_dqds);
    
    //We will used TSpectrum to get a list of peaks, problem is we have to put truncatedQ/dqds into TH1,
    //then make new TSpectrum, then run peakfinding then find the max "peaks" in each...
    std::swap(cluster._t_mean_v,truncated_mean);
    std::swap(cluster._t_dqds_v,truncated_dqds);
    
    if(candidate_loc == kINVALID_SIZE || dqdscandidate_loc == kINVALID_SIZE)
      return false;
    
    //20 is hardcoded
    if(abs(dqdscandidate_loc - candidate_loc) > 20)
      return false;
    
    auto window_size = 20;
    auto right = cluster._ordered_pts.size() - 1 - candidate_loc;
    auto left  = candidate_loc;
    
    
    bool right_is_smaller;
    int  iMin = 0;
    int  iMax = 0;
    
    if(right < left)
      right_is_smaller = true;
    else
      right_is_smaller = false;
    
    if(right > window_size) right = window_size;
    if(left  > window_size) left  = window_size;
    
    if(right_is_smaller) {
      iMax = right + candidate_loc;
      iMin = candidate_loc - window_size;
      if(iMin < 0) iMin = 0;
      
    } else {
      iMax = candidate_loc + window_size;
      iMin = candidate_loc - left;
      
      if(iMax >= cluster._ordered_pts.size()) iMax = cluster._ordered_pts.size() - 1;
    }
    
    auto k   = 0.0;
    auto idx = 0;

    for(int w = iMin; w <= iMax; ++w) {
      auto c = cluster.GetHits()[cluster._ordered_pts[w]]._q;
      if(c > k) { k = c; idx = w; }
    }
    
    
    cluster._boundary = cluster._ordered_pts[idx];
    return true;
  }
  
  
  unsigned int TSpectrumBoundary::nCk( unsigned int n, unsigned int k )
  {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
  
    int result = n;
    for( int i = 2; i <= k; ++i ) {
      result *= (n-i+1);
      result /= i;
    }
    return result;
  }
  

  double TSpectrumBoundary::coeff(double k, double N) {
    auto m = (N - 3.0)/2.0;
    return 1.0/pow(2,2*m+1) * (nCk(2*m,m-k+1) - nCk(2*m,m-k-1));
  }

  std::vector<double> TSpectrumBoundary::calc_smooth_derive(const std::vector<double>& _dist,
							     const std::vector<double>& tmeans, 
							     const int s)
  {
    std::vector<double> tdqds;
    tdqds.reserve(tmeans.size());
    
    if(!tmeans.size()) return tdqds;

    for(int o = 0; o < s; ++o) tdqds.push_back(0.0);
    
    //do smooth differentiation
    for(int i = s; i < tmeans.size() - s + 1; ++i) {
      std::vector<double> f(tmeans.begin() + i - s, tmeans.begin() + i + s);
      std::vector<double> x(_dist.begin() + i - s , _dist.begin() + i + s );
      tdqds.push_back(do_smooth_derive(f,x,2*s+1));
    }
    for(int o = 0; o < s - 1; ++o) tdqds.push_back(0.0);
    
    return tdqds;

  }
  
  
  double TSpectrumBoundary::do_smooth_derive(const std::vector<double>& f,
					      const std::vector<double>& x,
					      int N) 
  {
  
    // N should def be odd.
    auto M   = int{(N - 1)/2};
    auto tot = double{0.0};
  
    for(int k = 0; k < M; ++k)
      tot += coeff(k+1,N) * (f[k+M] - f[M - 1 - k])/(x[k + M] - x[M - 1 - k]) * 2 * (k+1);
  
    return tot;
  
  }
  
  std::vector<double> TSpectrumBoundary::calc_smooth_mean(const MichelCluster& cluster,
							   const double _n_window_size,
							   const int window_cutoff,
							   const double p_above) 
  {
    
    std::vector<double> charge;
    std::vector<double> truncatedQ;
    charge.reserve(cluster._ordered_pts.size());

    for(const auto& o : cluster._ordered_pts)
      charge.push_back(cluster.GetHits()[o]._q);

    for(auto window : get_windows(charge,_n_window_size) ) {
      if(window.size() > window_cutoff) 
	cut(window,p_above);
      truncatedQ.push_back(calc_mean(window));
    }
    return truncatedQ;
  }
  
  template<typename S>
  S TSpectrumBoundary::calc_mean(const std::vector<S>& data)
  {
    auto sum = S{0.0};
    for(const auto& d : data) sum += d;
    return sum / ( (S) data.size() ); //hopefully this isn't int :)
  }
  
  
  template<typename W>
  void TSpectrumBoundary::cut(std::vector<W>& data, double frac) 
  {
    
    auto size   = data.size();
    int to_stay = floor(frac*size);
  
    //sort the array based on charge
    std::sort(data.begin(),data.end(),
	      [](const W& a, const W& b) -> bool
	      {
		return a < b;	      
	      });
    
    //if(above) 
    data.erase(data.begin() + to_stay, data.end());
    //else 
    //data.erase(data.begin() , data.begin()+to_stay);
  }
  
  

  
  template<typename T>
  std::vector<std::vector<T> > TSpectrumBoundary::get_windows(const std::vector<T>& the_thing,
							       const int window_size)
  {
    
    std::vector<std::vector<T> > data;
    
    auto w = window_size + 2;
    w = (unsigned int)((w - 1)/2);
    auto num = the_thing.size();
    
    data.reserve(num);
    
    for(int i = 1; i <= num; ++i) {
      std::vector<T> inner;
      inner.reserve(20);
      if(i < w) {
	for(int j = 0; j < 2 * (i%w) - 1; ++j)
	  inner.push_back(the_thing[j]);
      }else if (i > num - w + 1){
	for(int j = num - 2*((num - i)%w)-1 ; j < num; ++j)
	  inner.push_back(the_thing[j]);
      }else{
	for(int j = i - w; j < i + w - 1; ++j)
	  inner.push_back(the_thing[j]);
      }
      data.emplace_back(inner);
    }

    return data;
  
  }
  
  size_t TSpectrumBoundary::find_max(const std::vector<double>& data) {
    
    auto the_max = double{0.0};
    size_t idx = kINVALID_SIZE;
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] > the_max) {
	the_max = data[i]; idx = i;
      }
    }
    
    return idx;
  }
  
  size_t TSpectrumBoundary::find_min(const std::vector<double>& data) {
  
    //get lowest dqds
    auto the_min = double{0.0};
    size_t idx= kINVALID_SIZE;
  
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] < the_min) {
	the_min = data[i]; idx = i;
      }
    }

    return idx;

  }

  size_t TSpectrumBoundary::get_tspectrum_max(const std::vector<double>& X, 
					      const std::vector<double>& Y) {
    
    double x[X.size() + 1];

    for(int i = 0; i <= X.size(); ++i) {
      if(i == X.size())x[i] = X.at(i-1) + 0.3;
      else             x[i] = X.at(i);

      
    }


    TSpectrum ts;                                                                                                        
    TH1D th("charge",";;",X.size(),x);    
    
    for(int i = 0; i < Y.size(); ++i) 
      th.SetBinContent(i+1,Y[i]); 
    
    ts.Search(&th,2,"nodraw");

    TH1D  *th_pt     = &th;
    TList *functions = th_pt->GetListOfFunctions(); 
    TPolyMarker* pm  = (TPolyMarker*)functions->FindObject("TPolyMarker"); 
    
    if(!pm)
      return kINVALID_SIZE;    

    std::vector<double> peakYs;
    peakYs.reserve(pm->GetN());
    
    for(int i = 0; i < pm->GetN(); ++i) {
      peakYs.push_back(pm->GetY()[i]);
    }
      
    auto maxY    = find_max(peakYs);
    
    if(maxY == kINVALID_SIZE)
      return kINVALID_SIZE;

    auto counter = int{-1};
    
    bool seen = false;
    for(const auto& y : Y) {
      counter++;
      if(y == peakYs.at(maxY)){
	seen = true;
	break;
      }
    }
    
    if(!seen)
      return kINVALID_SIZE;
    
    return (size_t)counter;
  }
  
  
}
#endif
