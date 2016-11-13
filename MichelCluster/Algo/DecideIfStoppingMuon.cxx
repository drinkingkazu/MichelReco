#ifndef DECIDEIFSTOPPINGMUON_CXX
#define DECIDEIFSTOPPINGMUON_CXX

#include "DecideIfStoppingMuon.h"
#include <sstream>
namespace michel{

  DecideIfStoppingMuon::DecideIfStoppingMuon()
  {
    _name    = "DecideIfStoppingMuon";
    _chi_min = 0.9;
    _frac_min_hits = 0.7;
    _hit_radius = 10;
    _max_dist   = 3;
    _min_bad_hits = 10;
  }

  bool DecideIfStoppingMuon::ProcessCluster(MichelCluster& cluster,
					    const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;

    // get the slope vector for the cluster
    auto const& slope_v = cluster._dirs_v;
    // get the chi vector
    auto const& chi_v = cluster._chi2_v;
    // get the hit vector
    auto const& hit_v = cluster._hits;
    // get the boundary point
    auto const& boundary = cluster._boundary;
    // is the michel forward (from boundary to end) or back
    auto const& forward = cluster._forward;
    // michel's start point
    auto const& start = cluster._hits[boundary];

    // how this algorithm works:
    // the idea is to take the straight track of a the muon-portion
    // of the cluster and get an average slope for the muon-like line
    // we anchor this slope to the average (w,t) point from the 
    // muon and thus get a line-equation for the muon-line
    // t = s * w + b
    // we then extend this line forward towards the michel direction
    // (using the point found to be the michel start point as our guide)
    // and look for ANY hits in the event that lie close to the michel
    // start point.
    // if these hits are well aligned with the muon-line (i.e. their
    // perpendicular distance to the line is small) then we tag these
    // as hits that actually belong to the muon.
    // if too many such hits exist, then we assume the muon actually
    // continues after the so called michel and therefore this is a
    // bad michel
    // this aims to resolve:
    // 1) deltas tagged as michels
    // 2) muons continuing after a section of broken wires


    // use the high-chi hits (value close to 1)
    // from the muon-section of the cluster
    // to find a direction for the muon itself
    // for these hits calculate an average slope
    // of the muon
    // project forward and check if there are many
    // other hits that fall along this slope
    // if so, probably we did not identify
    // correctly the end of the muon -> do not
    // reconstruct this as michel
    
    // average slope for the straight muon section:
    double slope = 0;
    int count= 0;
    // for those points that are "straight"
    // get the average x and z to identify
    // an anchor for the slope
    double w_avg = 0;
    double t_avg = 0;
    if (forward){
      for (size_t i=0; i < boundary; i++){
	if (fabs(chi_v[i]) > _chi_min){
	  slope += slope_v[i];
	  count += 1;
	  w_avg += hit_v[i]._w;
	  t_avg += hit_v[i]._t;
	}
      }
    }
    else{
      for (size_t i=boundary; i < chi_v.size(); i++){
	if (fabs(chi_v[i]) > _chi_min){
	  slope += slope_v[i];
	  count += 1;
	  w_avg += hit_v[i]._w;
	  t_avg += hit_v[i]._t;
	}
      }
    }
    
    // require that we have used at least 70% of the hits in the cluster
    double frac_used = (double)count / (double)hit_v.size();
    if(_verbosity <= msg::kINFO) {
      std::stringstream ss;
      ss << "Frac. of hits used to get slope: " << frac_used;
      Print(msg::kINFO,this->Name(),ss.str());
    }
    if (frac_used < _frac_min_hits)
	return false;

    slope /= count;
    w_avg /= count;
    t_avg /= count;
    // slope now is the average direction of the muon
    
    // now loop through all hits in the event
    // for those within some distance of the michel's
    // start point, if they fall close to the line
    // defined by y = s *x + b
    // where s is the slope, and b the michel start
    // coordinates, then count them as bad hits
    // potentially coming from the continuation of
    // the track.
    // if there are too many such hits -> BAD!


    // find anchor for slope (i.e. b in equation for line)
    // by using the "average muon point"
    double b = t_avg - slope * w_avg;
    
    if (_verbosity <=  msg::kINFO) {
      std::stringstream ss;
      ss << "\tDecideIfStoppingMuon Algo"<< std::endl
	 << "\tStart Point: [" << start._w << ", " << start._t << "]"<< std::endl
	 << "\tslope: " << slope<< std::endl
	 << "\thits used to calculate slope: " << count<< std::endl
	 << "\tline: y = s * x + b   ...  b = " << b ;
      Print(msg::kINFO,this->Name(),ss.str());
    }
    
    
    // hits identified as from continuation of muon:

    int nbad = 0;
    // loop through all hits
    for (auto const& h : hits){
      // if hit plane is not 2 -> skip
      if (h._pl != 2) continue;
      // distance to michel start
      double dMichel = h.SqDist(start);
      if ( dMichel < (_hit_radius*_hit_radius)){
	// if the distance to the michel is less then the distance to the center of the muon
	// then ok. Otherwise don't count the point because it is "behind" the michel
	//double dMuon = (h._w-w_avg)*(h._w-w_avg)+(h._t-t_avg)*(h._t-t_avg);
	//if (dMuon < dMichel) continue;
	// is this point in this cluster?
	bool use = true;
	for (auto const& hclus : cluster._hits){
	  if (h._id == hclus._id){
	    use = false;
	    break;
	  }
	}
	if (use){
	  // distance of this point to the line
	  double d = fabs(h._t - slope * h._w - b) / fabs(slope);
	  if (d < _max_dist) nbad += 1;
	}
      }
    }

    if (nbad > _min_bad_hits){
      if(_verbosity <= msg::kINFO) {
	std::stringstream ss;
	ss << "Bad hits aligned with muon: " << nbad;
	Print(msg::kINFO,this->Name(),ss.str());
      }
      return false;
    }

    return true;
  }

}

#endif
