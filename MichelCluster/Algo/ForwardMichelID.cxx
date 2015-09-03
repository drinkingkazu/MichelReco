#ifndef MICHELCLUSTER_FORWARDMICHELID_CXX
#define MICHELCLUSTER_FORWARDMICHELID_CXX

#include "ForwardMichelID.h"
#include "Fmwk/MichelException.h"

namespace michel {

  void ForwardMichelID::EventReset()
  {}
  
  bool ForwardMichelID::ProcessCluster(MichelCluster& cluster,
				       const std::vector<HitPt>& hits)
  {
    
    if(cluster._boundary == kINVALID_SIZE) { return false; }
    //no check on cluster size here...
    
    //Hardcoded for now :)
    double n_cutoff = 2;
    double c_cutoff = 1.15;
    int    w_cutoff = 10;
    int    idx = -1;
    size_t the_michel_start; //index in ordered points that is the michel
    
    //get index in ordered_pts that is the boundary point
    
    for(const auto& ordered_hitIDX : cluster._ordered_pts) {
      idx++;
      if(cluster._boundary == ordered_hitIDX) break;
    }
    
    //Determines which point is the michel start    
    if(determine_forward(cluster,n_cutoff,c_cutoff,w_cutoff,cluster._forward)) {
      if(cluster._forward) {
	if(idx >= cluster._ordered_pts.size() - 1) {
	  the_michel_start = idx;
	}
	else {
	  the_michel_start = idx + 1;
	}
      } 
      //forward == false
      else {
	if(idx != 0) {
	  the_michel_start = idx - 1;
	}
	else {
	  the_michel_start = idx;
	}
      } 
    }
    
    //we could not determine forward so we return default contructor. This is a "bad" michel.   
    else { return false; }
    
    auto& electron = cluster._michel;
    std::vector<size_t> ordered_pts_idx;

    for( size_t i = 0 ; i < cluster._ordered_pts.size(); ++i) {
      if(cluster._forward)  {
	if(i >= the_michel_start) {
	  electron.push_back(cluster._hits[cluster._ordered_pts[i]]);
	  ordered_pts_idx.push_back(i);
	}
      }
      else {
	if(i <= the_michel_start)  {
	  electron.push_back(cluster._hits[cluster._ordered_pts[i]]);
	  ordered_pts_idx.push_back(i);
	}
      }
    }

    // if the michel has more than _maxHits points -> do not return a michel...it is probably garbage
    if ( (electron.size() > _maxHits) and (_maxHits != 0) )
      return false;

    //Do same thing for muon (in future)
    auto length = determine_length(cluster,ordered_pts_idx); //true radius with no minimum
    
    electron._length = length;
    electron._start  = cluster._hits[cluster._ordered_pts[the_michel_start]];
    
    //loop over the ordered points, add hits to the electron that are NOT in orderedpts
    //BUT append new hits found in the vicinity of the vertex
    //this now goes into "reclustering"
    return true;
  }

  bool ForwardMichelID::determine_forward(const MichelCluster& cluster,
					  const double n_cutoff,
					  const double c_cutoff,
					  const double w_cutoff,
					  bool&  forward) 
  {
    
    
    double lower_idx_Q  = 0;
    double higher_idx_Q = 0;
    
    unsigned int lower_idx  = 0;
    unsigned int higher_idx = 0;
    
    unsigned int N_ordered_pts = cluster._ordered_pts.size();

    // Loop over the ordered cluster points, and add up the number of points left
    // and right of the boundary. Also add up the total Q in the same manner.

    for (HitIdx_t i = 0; i < N_ordered_pts; ++i){
      if (i < cluster._boundary) {
	lower_idx++;
	lower_idx_Q += cluster._hits[cluster._ordered_pts[i]]._q;
      }
      else if (i > cluster._boundary){
	higher_idx++;
	higher_idx_Q += cluster._hits[cluster._ordered_pts[i]]._q;
      }
    }
    
    if(higher_idx == 0) { forward = true;  return true; }
    if(lower_idx  == 0) { forward = false; return true; }
    if(n_cutoff   == 0) { throw MichelException();       }
    if(c_cutoff   == 0) { throw MichelException();       }
    
    //First, check if the ratio of # of points on left and right is above specified cutoff..
    if( lower_idx/higher_idx > n_cutoff || lower_idx/higher_idx < 1/n_cutoff) {
      
      if (lower_idx_Q > higher_idx_Q)
	forward = true;
      else
	forward = false; 
    }
    
    //Second, check if the ratio of the charges on the let and right is above specified cutoff...
    else if (lower_idx_Q/higher_idx_Q > c_cutoff || lower_idx_Q/higher_idx_Q < 1/c_cutoff) {
      
      if (lower_idx_Q > higher_idx_Q) 
	forward = true;
      else
	forward = false; 
    }
    
    //Third, look in window around boundary points, see if one side has more charge or 
    else if(lower_idx > w_cutoff && higher_idx > w_cutoff) {
      
      double lower_idx_Q_window  = 0.0;
      double higher_idx_Q_window = 0.0;
      
      for (size_t i = cluster._boundary - w_cutoff; i < cluster._boundary + w_cutoff; ++i){
    	if (i < cluster._boundary)
    	  lower_idx_Q_window  += cluster._hits[cluster._ordered_pts[i]]._q;
    	else if (i > cluster._boundary) 
    	  higher_idx_Q_window += cluster._hits[cluster._ordered_pts[i]]._q;
      }
      
      if (lower_idx_Q_window > higher_idx_Q_window)
	forward = true;
      else
	forward = false; 
    }
    
    else { //no good, need new ALGO?
      forward = false;
      return false;  
    }
    
    return true;
  }
  
  double ForwardMichelID::determine_length(const MichelCluster& c,
					   const std::vector<size_t>& ordered_pts_idx) 
  {
    double length = 0;
    
    for(const auto& o : ordered_pts_idx) 
      if(o < c._ds_v.size())
	length += c._ds_v[o];
    
    return length;
  }
  
}
#endif
