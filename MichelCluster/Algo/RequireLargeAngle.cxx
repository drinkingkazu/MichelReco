#ifndef REQUIRELARGEANGLE_CXX
#define REQUIRELARGEANGLE_CXX

#include "RequireLargeAngle.h"
#include "Fmwk/ClusterVectorCalculator.h"
#include <sstream>
#include <cmath>

namespace michel{

  RequireLargeAngle::RequireLargeAngle()
  {
    _name    = "RequireLargeAngle";
    _min_angle = (20.*3.14/180.);
    _min_straight_michel_hits = 0;
    _muon_length_used = 20;
  }

  bool RequireLargeAngle::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    /// call instance of "ClusterVectorCalculator"
    /// this class has a bunch of utility functions
    /// to calculate stuff based on the vecor of
    /// hits for a cluster
    ClusterVectorCalculator _clusterCalc;


    // This algorithm tries to find a slope
    // for the muon and one for the michel
    // that indicate the general direction
    // in a 2D plane for these two "tracks"
    // we require that these tracks are 
    // at a large enough angle one w.r.t.
    // the other.
    // otherwise, we remove the event

    auto const& michel_hits = cluster._michel;

    auto const& michel_slope = _clusterCalc.calc_slope(michel_hits,5);
    auto const& michel_chi   = _clusterCalc.calc_covariance(michel_hits,5);

    if(_verbosity <= msg::kINFO) {
      std::cout << std::endl;
      std::cout << "Michel Slope vector: [";
      for (auto& s : michel_slope){
	printf("%.02f",s);
	std::cout << ", ";
      }
      std::cout << "]" << std::endl;
      std::cout << "Michel Chi vector: [";
      for (auto& c : michel_chi){
	printf("%.02f",c);
	std::cout << ", ";
      }
      std::cout << "]" << std::endl;
    }

    // calculate the average slope for the michel
    int count = 0;
    double avg_slope  = 0.;
    for (size_t n=0; n < michel_chi.size(); n++){
      if (fabs(michel_chi[n]) > 0.7){
	avg_slope += michel_slope[n];
	count += 1;
      }
    }// for all points in michel cluster
    if (count == 0){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "++++++++ We found 0 hits from which to calculate the slope of the michel -> rejecting this event...";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }

    // require the michel to have a minimum number of "straight hits"
    if (count < _min_straight_michel_hits){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "++++++++ number of michel hits with chi < " << 0.5 << " is too small [" << count << "]. Removing event...";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }
    
    double slope_michel = avg_slope/((double)count);

    // David C Aug 9th 2016

    auto const& michel_dir = cluster._michel._dir;

    slope_michel = michel_dir._t / michel_dir._w;

    // same for the muon -> calculate an average slope

    // get relevant quantities
    // get the slope vector for the cluster
    auto const& slope_v = cluster._dirs_v;
    // get the chi vector
    auto const& chi_v = cluster._chi2_v;
    // get the dS vector
    auto const& dist_v = cluster._s_v;
    // get the boundary point
    auto const& boundary = (int)cluster._boundary;
    // is the michel forward (from boundary to end) or back
    auto const& forward = cluster._forward;

    // average slope for the straight muon section:
    if (forward){
      for (size_t i=0; i < boundary; i++){
	// only consider points w/ distance to boundary < _muon_length_used
	if ( fabs( dist_v[i] - dist_v[boundary] ) > _muon_length_used ) continue;
	if (fabs(chi_v[i]) < 0.7) continue;
	avg_slope += slope_v[i];
	count += 1;
      }
    }
    else{
      for (size_t i=boundary; i < chi_v.size(); i++){
	if ( fabs( dist_v[i] - dist_v[boundary] ) > _muon_length_used ) continue;
	if (fabs(chi_v[i]) > 0.7) continue;
	avg_slope += slope_v[i];
	count += 1;
      }
    }

    if (count == 0){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "++++++++ We found 0 hits from which to calculate the slope of the muon -> rejecting this event...";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }

    double slope_muon = avg_slope/((double)count);

    // finally, get the angle between the two slopes
    
    double arctan = (slope_michel - slope_muon) / ( 1 + slope_michel * slope_muon );

    double angle = fabs(atan(arctan));

    if(_verbosity <= msg::kINFO){
      std::cout << "Muon slope   : " << slope_muon << std::endl;
      std::cout << "Michel slope : " << slope_michel << std::endl;
      std::cout << "Angle is: " << angle << std::endl;
    }

    if (angle < _min_angle){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "++++++++ angle b/w muon and michel is too small. Found ot be " << angle << ". rejecting event...";
	Print(msg::kINFO,__FUNCTION__,ss.str());
      }
      return false;
    }

    return true;
  }

}

#endif
