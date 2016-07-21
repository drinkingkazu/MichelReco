#ifndef CLUSTERPHOTONS_CXX
#define CLUSTERPHOTONS_CXX

#include "ClusterPhotons.h"

namespace michel{

  ClusterPhotons::ClusterPhotons()
  {
    _name    = "ClusterPhotons";
    _d_max   = 3;
  }

  bool ClusterPhotons::ProcessCluster(MichelCluster& cluster,
				      const std::vector<HitPt>& hits)
  {

    // grab michel hits
    auto const& michel_hits = cluster._michel;

    // clear the map that links hit index to cluster they have been saved to
    _hit_cluster_map.clear();

    // clear this michel's cluster of photons
    cluster._michel._photon_clus_v.clear();

    // loop trhough all hits in michel
    for (size_t i=0; i < michel_hits.size(); i++){

      auto const& pt1 = michel_hits[i];

      for (size_t j=i+1; j < michel_hits.size(); j++){

	auto const& pt2 = michel_hits[j];

	// if close enough
	if (pt1.SqDist(pt2) < _d_max) {

	  // is hit i recorded in the map?
	  if ( _hit_cluster_map.find( i ) != _hit_cluster_map.end() ){
	    cluster._michel._photon_clus_v[ _hit_cluster_map[ i ] ].push_back( j );
	    _hit_cluster_map[j] = _hit_cluster_map[i];
	  }
	  // is hit j recorded in the map?
	  else if ( _hit_cluster_map.find( j ) != _hit_cluster_map.end() ){
	      cluster._michel._photon_clus_v[ _hit_cluster_map[ j ] ].push_back( i );
	      _hit_cluster_map[i] = _hit_cluster_map[j];
	  }
	  else{
	    std::vector<size_t> new_clus = {i,j};
	    cluster._michel._photon_clus_v.push_back( new_clus );
	    _hit_cluster_map[i] = cluster._michel._photon_clus_v.size() - 1;
	    _hit_cluster_map[j] = cluster._michel._photon_clus_v.size() - 1;
	  }
	  
	}// if points are close
	
      }// 2nd loop throguh hits
    }// 1st loop through hits
    
    
    if(_verbosity <= msg::kINFO) {
      std::stringstream ss;
      ss << "\n\t\tFound " << cluster._michel._photon_clus_v.size()
	 << " sub-clusters in Michel."<< std::endl;
      Print(msg::kINFO,__FUNCTION__,ss.str());
    }
    
    return true;
  }

}

#endif
