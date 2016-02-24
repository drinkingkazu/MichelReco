#ifndef MICHELCLUSTER_CLUSTERVECTORCALCULATOR_CXX
#define MICHELCLUSTER_CLUSTERVECTORCALCULATOR_CXX

#include "ClusterVectorCalculator.h"
#include "Fmwk/MichelException.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <sstream>

namespace michel {
  
  unsigned int ClusterVectorCalculator::nCk( unsigned int n, unsigned int k ) const
  {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
  
    int result = n;
    for( unsigned int i = 2; i <= k; ++i ) {
      result *= (n-i+1);
      result /= i;
    }
    return result;
  }
  

  double ClusterVectorCalculator::coeff(double k, double N)  const
  {
    auto m = (N - 3.0)/2.0;
    return 1.0/pow(2,2*m+1) * (nCk(2*m,m-k+1) - nCk(2*m,m-k-1));
  }

  std::vector<double> ClusterVectorCalculator::calc_smooth_derive(const std::vector<double>& _dist,
							     const std::vector<double>& tmeans, 
							     const int s) const
  {
    std::vector<double> tdqds;
    tdqds.reserve(tmeans.size());
    
    if(!tmeans.size()) return tdqds;

    for(int o = 0; o < s; ++o) tdqds.push_back(0.0);
    
    //do smooth differentiation
    for(int i = s; i < (int)tmeans.size() - s + 1; ++i) {
      std::vector<double> f(tmeans.begin() + i - s, tmeans.begin() + i + s);
      std::vector<double> x(_dist.begin() + i - s , _dist.begin() + i + s );
      tdqds.push_back(do_smooth_derive(f,x,2*s+1));
    }
    for(int o = 0; o < s - 1; ++o) tdqds.push_back(0.0);
    
    return tdqds;

  }
  
  
  double ClusterVectorCalculator::do_smooth_derive(const std::vector<double>& f,
					      const std::vector<double>& x,
					      int N) const
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
								const double p_above) const
  {
    
    std::vector<double> charge;
    std::vector<double> truncatedQ;
    charge.reserve(cluster._ordered_pts.size());

    for(const auto& o : cluster._ordered_pts)
      charge.push_back(cluster._hits[o]._q);

    for(auto window : get_windows(charge,_n_window_size) ) {
      if(window.size() > (size_t)window_cutoff) 
	cut(window,p_above);
      truncatedQ.push_back(calc_mean(window));
    }
    return truncatedQ;
  }
  
  template<typename S>
  S ClusterVectorCalculator::calc_mean(const std::vector<S>& data) const
  {
    if(!data.size())
      Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to calc_mean");

    auto sum = S{0.0};
    for(const auto& d : data) sum += d;
    return sum / ( (S) data.size() ); //hopefully this isn't int :)
  }
  
  
  template<typename W>
  void ClusterVectorCalculator::cut(std::vector<W>& data, double frac) const
  {
    
    auto size   = data.size();
    // calcualte number of elements to be kept
    int to_stay = floor(frac*size);
  
    // sort the array based on charge
    // so that high-charge hits are removed
    std::sort(data.begin(),data.end(),
	      [](const W& a, const W& b) -> bool
	      {
		return a < b;	      
	      });

    // erase all elements after the last one to be kept
    data.erase(data.begin(), data.begin() + to_stay);
    data.erase(data.end() - to_stay, data.end());
    //data.erase(data.begin() + to_stay, data.end());
  }
  
  

  
  template<typename T>
  std::vector<std::vector<T> > ClusterVectorCalculator::get_windows(const std::vector<T>& the_thing,
							       const int window_size) const
  {

    // given a vector of values return a vector of the same length
    // with each element being a vector of the values of the local neighbors
    // of the element at position i in the original vector
    // input  : [0,1,2,3,4,5,6,...,...,N-3,N-2,N-1] (input vector of size N)
    // output  (assuming a value of 'w' below == 3):
    // 0th element: [0]
    // 1st element: [0,1,2]
    // 2nd element: [0,1,2,3,4]
    // jth element: [j-w,j-w+1,..,j+w-2,j+w-1]
    
    std::vector<std::vector<T> > data;
    
    auto w = window_size + 2;
    w = (unsigned int)((w - 1)/2);
    auto num = the_thing.size();
    
    data.reserve(num);
    
    for(int i = 1; i <= num; ++i) {
      std::vector<T> inner;
      inner.reserve(20);
      // if we are at the beginning of the vector (and risk accessing -1 elements)
      if(i < w)
	{
	  for(int j = 0; j < 2 * (i%w) - 1; ++j)
	    inner.push_back(the_thing[j]);
	}
      // if we are at the end of the vector (and risk going past it)
      else if (i > num - w + 1)
	{
	  for(int j = num - 2*((num - i)%w)-1 ; j < num; ++j)
	    inner.push_back(the_thing[j]);
	}
      // if we are in the middle of the waveform
      else
	{
	  for(int j = i - w; j < i + w - 1; ++j)
	    inner.push_back(the_thing[j]);
	}
      data.emplace_back(inner);
    }

    return data;
  
  }
  
  size_t ClusterVectorCalculator::find_max(const std::vector<double>& data) const
  {
    
    auto the_max = double{0.0};
    size_t idx = kINVALID_SIZE;
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] > the_max) {
	the_max = data[i]; idx = i;
      }
    }
    
    return idx;
  }
  
  size_t ClusterVectorCalculator::find_min(const std::vector<double>& data) const
  {
  
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
							  const int _n_window_size) const
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
							   const int _n_window_size) const
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
				  const std::vector<double>& data2) const
  {
    if(data1.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to cov");
    if(data2.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to cov");

    double result = 0.0;
    auto   mean1  = mean(data1);
    auto   mean2  = mean(data2);
    
    for(size_t i = 0; i < data1.size(); ++i)
      result += (data1[i] - mean1)*(data2[i] - mean2);
    
    return result/((double)data1.size());
      
  }
  
  double ClusterVectorCalculator::stdev(const std::vector<double>& data) const
  {
    if(data.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to stdev");

    double result = 0.0;
    auto    avg   = mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));
  }
  
  double ClusterVectorCalculator::mean(const std::vector<double>& data) const
  {
    if(data.size() == 0) Print(msg::kEXCEPTION,__FUNCTION__,"You have me nill to mean");
	
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
        
    return (result / ((double)data.size()));
  }


  std::pair<double,double> ClusterVectorCalculator::GetLinearFit(const std::vector<michel::HitPt>& pts) const
  {

    // mean of x*y
    double m_xy = 0;
    // mean of x
    double m_x  = 0;
    // mean of y
    double m_y  = 0;
    // mean of x^2
    double m_xx = 0;
    // mean of y^2
    double m_yy = 0;
    
    for (size_t n=0; n < pts.size(); n++){

      double xi = pts[n]._w;
      double yi = pts[n]._t;

      m_xy += xi*yi;
      m_x  += xi;
      m_y  += yi;
      m_xx += xi*xi;
      m_yy += yi*yi;

    }

    double entries = (double)pts.size();

    m_xy /= entries;
    m_x  /= entries;
    m_y  /= entries;
    m_xx /= entries;
    m_yy /= entries;

    double slope = (m_xy-m_x*m_y)/(m_xx-m_x*m_x);
    double intercept = m_y - slope * m_x;

    return std::pair<double,double>(slope,intercept);
  }


  std::pair<double,double> ClusterVectorCalculator::GetLinearFit(const std::vector<double>& x, 
								 const std::vector<double>& y)
  {

    // make sure both vectors are the same size!
    if (x.size() != y.size()){
      std::cerr << "\033[93m[ERROR]\033[00m dS and dQ vectors do not have the same length"
		<< std::endl;
      throw michel::MichelException();
    }

    // mean of x*y
    double m_xy = 0;
    // mean of x
    double m_x  = 0;
    // mean of y
    double m_y  = 0;
    // mean of x^2
    double m_xx = 0;
    // mean of y^2
    double m_yy = 0;
    
    for (size_t n=0; n < x.size(); n++){

      double xi = x[n];
      double yi = y[n];

      m_xy += xi*yi;
      m_x  += xi;
      m_y  += yi;
      m_xx += xi*xi;
      m_yy += yi*yi;

    }

    double entries = (double)x.size();

    m_xy /= entries;
    m_x  /= entries;
    m_y  /= entries;
    m_xx /= entries;
    m_yy /= entries;

    double slope = (m_xy-m_x*m_y)/(m_xx-m_x*m_x);
    double intercept = m_y - slope * m_x;

    return std::pair<double,double>(slope,intercept);
  }


  double ClusterVectorCalculator::GetPerpendicularDistance(const michel::HitPt& h,
							   const double& slope,
							   const double& intercept) const
  {

    if (slope == 0)
      return fabs ( h._t - intercept );

    double d_perp = fabs ( h._t - slope * h._w - intercept ) / fabs(slope);

    return d_perp;
  }

  // get the median value for a list of doubles
  double ClusterVectorCalculator::GetMedian(std::vector<double>& v) const
  {

    std::nth_element(v.begin(),v.begin() + v.size()/2, v.end());

    return v[v.size()/2];
  }

  size_t ClusterVectorCalculator::GetMIPendPos(const std::vector<double>& v,
					       const size_t& max,
					       const double distAsked) const
  {

    double distToPeak = v[max];
    
    for (size_t i=0; i < max; i++){
      
      double dist = distToPeak-v[max-i];
      if (dist > distAsked)
	return max-i;
    }

    //std::cout << "could not find a point " << distAsked
    //      << " away from the bragg peak..." << std::endl;

    return michel::kINVALID_SIZE;
  }


  std::vector<size_t> ClusterVectorCalculator::GetMIPindices(const std::vector<double>& dQ,
							     const double& median,
							     const double& rms,
							     const double& alpha) const
  {

    std::vector<size_t> indices;

    for (size_t n=0; n < dQ.size(); n++){
      if ( (dQ[n] < median+rms) and (dQ[n] > median-rms) )
	indices.push_back(n);
    }

    return indices;
  }


  std::vector<double> ClusterVectorCalculator::GetSubVector(const std::vector<double> v,
							    const std::vector<size_t> idx) const
  {

    std::vector<double> out;

    for (auto& i : idx){
      if (i >= v.size()){
	std::cerr << "\033[93m[ERROR]\033[00m trying to access element larger than vector size..."
		  << std::endl;
	throw michel::MichelException();
      }
      out.push_back(v[i]);
    }
    
    return out;
  }

  double ClusterVectorCalculator::GetRms(const std::vector<double>& MIPdS,
					 const std::vector<double>& MIPdQ,
					 const double& slope,
					 const double& intercept) const
  {

    double rms = 0.;
    
    for (size_t i=0; i < MIPdS.size(); i++)
      rms += (MIPdQ[i] - ( intercept + slope * MIPdS[i] )) * (MIPdQ[i] - ( intercept + slope * MIPdS[i] ));
    
    if (MIPdS.size() == 0)
      return 0;
  
    rms = sqrt(rms/(double)MIPdS.size());
    
    return rms;
  }


  std::pair<size_t,double> ClusterVectorCalculator::GetMaxIndex(const std::vector<double>& v) const
  {

    double max = 0;
    size_t idx = 0;
    for (size_t i=0; i < v.size(); i++)
      if (v[i] > max) { max = v[i]; idx = i; }

    return std::pair<size_t,double>(idx,max);
  }  


  double ClusterVectorCalculator::GetBraggArea(const std::vector<double>& dS,
					       const std::vector<double>& dQ,
					       const size_t& MIPendIdx,
					       const size_t& braggIdx,
					       const double& MIPm,
					       const double& MIPs) const
  {

    if (MIPendIdx >= braggIdx){
      std::cerr << "\033[93m[ERROR]\033[00m MIP end index comes after Bragg Peak index...something is wrong..."
		<< std::endl;
      throw michel::MichelException();
    }

    if (braggIdx > dS.size()){
      std::cerr << "\033[93m[ERROR]\033[00m Bragg peak index is larger than vector size...something is wrong..."
		<< std::endl;
      throw michel::MichelException();
    }

    if (dQ.size() != dS.size()){
      std::cerr << "\033[93m[ERROR]\033[00m dS and dQ vectors are not the same length...something is wrong..."
		<< std::endl;
      throw michel::MichelException();
    }

    double area = 0;
    for (size_t i =  MIPendIdx; i < braggIdx; i++)
      area += dQ[i] - ( MIPm + MIPs * dS[i] );

    return area;
  }


}
#endif
