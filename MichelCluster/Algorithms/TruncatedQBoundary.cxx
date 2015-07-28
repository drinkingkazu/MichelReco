#ifndef MICHELCLUSTER_TRUNCATEDQBOUNDARY_CXX
#define MICHELCLUSTER_TRUNCATEDQBOUNDARY_CXX

#include "TruncatedQBoundary.h"
#include "Fmwk/MichelException.h"
#include <cmath>

namespace michel {
  
  void TruncatedQBoundary::EventReset()
  {}
  
  HitIdx_t TruncatedQBoundary::Boundary(MichelCluster& cluster)
  { 
    
    std::vector<double> truncated_mean;
    std::vector<double> truncated_dqds;
    std::vector<double> covariance; 
    
    truncated_mean.reserve(cluster._ordered_pts.size());
    truncated_dqds.reserve(cluster._ordered_pts.size());
    covariance.reserve    (cluster._ordered_pts.size());
    
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
    covariance = calc_covariance(cluster._hits,11);
    
    //Lets play with truncated mean shaving...
    if(_verbosity <= msg::kINFO) {
      std::cout << "\n\t\tIn TruncatedQBoundary\n"
		<< "\tI have " << truncated_mean.size() << " truncated mean size\n"
		<< "\twith   " << truncated_dqds.size() << " derivative points.\n"
		<< "\tMy incoming cluster has " << cluster._hits.size() << " hits in it...\n";
    }
    
    //With this new information, calculate the boundary point between possible muon end and michel start
    
    auto candidate_loc     = find_max(truncated_mean);
    auto dqdscandidate_loc = find_min(truncated_dqds); 

    std::swap(cluster._t_mean_v,truncated_mean);
    std::swap(cluster._t_dqds_v,truncated_dqds);
    
    if((candidate_loc     >= cluster._hits.size()))
      return kINVALID_SIZE;
    
    if((dqdscandidate_loc >= cluster._hits.size()))
      return kINVALID_SIZE;
    
    if(abs(dqdscandidate_loc - candidate_loc) > 20)
      return kINVALID_SIZE;
    

    //20 is hardcoded
    auto window_size = 20;
    auto right = cluster._ordered_pts.size() - 1 - candidate_loc;
    auto left  = candidate_loc;
    
    int  iMin = 0;
    int  iMax = 0;
    
    
    if(right >= window_size) iMax  = window_size   + candidate_loc;
    if(left  >= window_size) iMin  = candidate_loc - window_size;

    if(right < window_size)  iMax  = cluster._hits.size() - 1;
    if(left  < window_size)  iMin  = 0;
    
    auto k   = 0.0;
    auto idx = 0;
    
    for(int w = iMin; w <= iMax; ++w) {
      auto c = cluster._hits[cluster._ordered_pts[w]]._q;
      if(c > k) { k = c; idx = w; }
    }
    
    std::swap(cluster._chi2_v,covariance);

    return cluster._ordered_pts[idx];
  }
  
  unsigned int TruncatedQBoundary::nCk( unsigned int n, unsigned int k )
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
  

  double TruncatedQBoundary::coeff(double k, double N) {
    auto m = (N - 3.0)/2.0;
    return 1.0/pow(2,2*m+1) * (nCk(2*m,m-k+1) - nCk(2*m,m-k-1));
  }

  std::vector<double> TruncatedQBoundary::calc_smooth_derive(const std::vector<double>& _dist,
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
  
  
  double TruncatedQBoundary::do_smooth_derive(const std::vector<double>& f,
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
  
  std::vector<double> TruncatedQBoundary::calc_smooth_mean(const MichelCluster& cluster,
							   const double _n_window_size,
							   const int window_cutoff,
							   const double p_above) 
  {
    
    std::vector<double> charge;
    std::vector<double> truncatedQ;
    charge.reserve(cluster._ordered_pts.size());

    for(const auto& o : cluster._ordered_pts)
      charge.push_back(cluster._hits[o]._q);

    for(auto window : get_windows(charge,_n_window_size) ) {
      if(window.size() > window_cutoff) 
	cut(window,p_above);
      truncatedQ.push_back(calc_mean(window));
    }
    return truncatedQ;
  }
  
  template<typename S>
  S TruncatedQBoundary::calc_mean(const std::vector<S>& data)
  {
    if(!data.size()){ std::cout << "You have me nill to calc_mean\n"; throw MichelException(); }
    auto sum = S{0.0};
    for(const auto& d : data) sum += d;
    return sum / ( (S) data.size() ); //hopefully this isn't int :)
  }
  
  
  template<typename W>
  void TruncatedQBoundary::cut(std::vector<W>& data, double frac) 
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
  std::vector<std::vector<T> > TruncatedQBoundary::get_windows(const std::vector<T>& the_thing,
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
  
  size_t TruncatedQBoundary::find_max(const std::vector<double>& data) {
    
    auto the_max = double{0.0};
    size_t idx = kINVALID_SIZE;
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] > the_max) {
	the_max = data[i]; idx = i;
      }
    }
    
    return idx;
  }
  
  size_t TruncatedQBoundary::find_min(const std::vector<double>& data) {
  
    //get lowest dqds
    auto the_min = double{999999.0};
    size_t idx= kINVALID_SIZE;
  
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] < the_min) {
	the_min = data[i]; idx = i;
      }
    }

    return idx;

  }

  std::vector<double>  TruncatedQBoundary::calc_covariance(const std::vector<::michel::HitPt>& hits, 
							   const int _n_window_size)
  {
    std::vector<double> R;
    R.reserve(hits.size());

    std::vector<double> X;
    std::vector<double> Y;
    X.reserve(_n_window_size);
    Y.reserve(_n_window_size);

    for(const auto& window : get_windows(hits,_n_window_size) ) {
      for(const auto& hit : window) {
	X.push_back(hit._w); Y.push_back(hit._t);
      }
      
      auto c  = cov(X,Y);
      auto sX = stdev(X);
      auto sY = stdev(Y);

      R.push_back(c / ( sX * sY));

      X.clear(); Y.clear();
    }    
    
    return R;
  }
  
  
  double TruncatedQBoundary::cov (const std::vector<double>& data1,
				  const std::vector<double>& data2)
  {
    if(!data1.size()){ std::cout << "You have me nill to cov\n"; throw MichelException(); }
    if(!data2.size()){ std::cout << "You have me nill to cov\n"; throw MichelException(); }

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(int i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
  }
  
  double TruncatedQBoundary::stdev(const std::vector<double>& data)
  {
    if(!data.size()){ std::cout << "You have me nill to stdev\n"; throw MichelException(); }

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }
  
  double TruncatedQBoundary::mean(const std::vector<double>& data)
  {
    if(!data.size()){ std::cout << "You have me nill to mean\n"; throw MichelException(); }
	
    double result = 0.0;
    for(const auto& d : data)
      result += d;

    return (result / ((double)data.size()));
  }


}
#endif
