#ifndef MICHELCLUSTER_CLUSTERVECTORCALCULATOR_CXX
#define MICHELCLUSTER_CLUSTERVECTORCALCULATOR_CXX

#include "ClusterVectorCalculator.h"
#include "Fmwk/MichelException.h"
#include <cmath>
#include <sstream>
namespace michel {
  
  unsigned int ClusterVectorCalculator::nCk( unsigned int n, unsigned int k )
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
  

  double ClusterVectorCalculator::coeff(double k, double N) {
    auto m = (N - 3.0)/2.0;
    return 1.0/pow(2,2*m+1) * (nCk(2*m,m-k+1) - nCk(2*m,m-k-1));
  }

  std::vector<double> ClusterVectorCalculator::calc_smooth_derive(const std::vector<double>& _dist,
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
  
  
  double ClusterVectorCalculator::do_smooth_derive(const std::vector<double>& f,
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
  
  std::vector<double> ClusterVectorCalculator::calc_smooth_mean(const MichelCluster& cluster,
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
  S ClusterVectorCalculator::calc_mean(const std::vector<S>& data)
  {
    if(!data.size())
      Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to calc_mean");

    auto sum = S{0.0};
    for(const auto& d : data) sum += d;
    return sum / ( (S) data.size() ); //hopefully this isn't int :)
  }
  
  
  template<typename W>
  void ClusterVectorCalculator::cut(std::vector<W>& data, double frac) 
  {
    
    auto size   = data.size();
    int to_stay = floor(frac*size);
  
    //sort the array based on charge
    std::sort(data.begin(),data.end(),
	      [](const W& a, const W& b) -> bool
	      {
		return a < b;	      
	      });
    
    data.erase(data.begin() + to_stay, data.end());
  }
  
  

  
  template<typename T>
  std::vector<std::vector<T> > ClusterVectorCalculator::get_windows(const std::vector<T>& the_thing,
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
  
  size_t ClusterVectorCalculator::find_max(const std::vector<double>& data) {
    
    auto the_max = double{0.0};
    size_t idx = kINVALID_SIZE;
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] > the_max) {
	the_max = data[i]; idx = i;
      }
    }
    
    return idx;
  }
  
  size_t ClusterVectorCalculator::find_min(const std::vector<double>& data) {
  
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

  std::vector<double> ClusterVectorCalculator::calc_slope(const std::vector<::michel::HitPt>& hits,
							  const int _n_window_size)
  {
    std::vector<double> S;
    S.reserve(hits.size());

    std::vector<double> X;
    std::vector<double> Y;
    X.reserve(_n_window_size);
    Y.reserve(_n_window_size);

    for(const auto& window : get_windows(hits,_n_window_size) ) { 
      for(const auto& hit : window) {
	X.push_back(hit._w); Y.push_back(hit._t);
      }
      
      //http://mathworld.wolfram.com/LeastSquaresFitting.html
      auto c   = cov(X,Y);
      auto sX  = stdev(X);
      if(sX == 0.0) {c = 0.0; sX = 1.0;} //might have a problem with floating point error here
      auto b   = c/(sX*sX) ;
      
      S.push_back(b);
      X.clear(); Y.clear();
    }
    S.at(0)            = S.at(1);
    S.at(S.size() - 1) = S.at(S.size() - 2);
    
    return S;
  }

  std::vector<double>  ClusterVectorCalculator::calc_covariance(const std::vector<::michel::HitPt>& hits, 
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
      auto r  = c/(sX * sY);

      if(_verbosity <= msg::kDEBUG) {
	std::stringstream ss;
	ss << "c: "  << c << std::endl
	   << "sX: " << sX <<  std::endl
	   << "sY: " << sY <<  std::endl
	   << "r: "  << r <<  std::endl;
	Print(msg::kDEBUG,__FUNCTION__,ss.str());
      }
      if(isnan(r)) r = 0.0; 
      R.push_back(r);
      

      // if(R.size() != 1 && R.size() != hits.size())
      // 	if(isnan(r)) Print(msg::kEXCEPTION,__FUNCTION__,"Covariance is nan not on edge");
      
      X.clear(); Y.clear();
    }    
    //first and last points will be nan. Lets set them equal to the points just above and below
    R.at(0)            = R.at(1);
    R.at(R.size() - 1) = R.at(R.size() - 2);
    
    return R;
  }
  
  
  double ClusterVectorCalculator::cov (const std::vector<double>& data1,
				  const std::vector<double>& data2)
  {
    if(data1.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to cov");
    if(data2.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to cov");

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(int i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
      
  }
  
  double ClusterVectorCalculator::stdev(const std::vector<double>& data)
  {
    if(data.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to stdev");

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }
  
  double ClusterVectorCalculator::mean(const std::vector<double>& data)
  {
    if(data.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to mean");
	
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }


}
#endif
