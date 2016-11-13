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

    if (hits.size() == 0) return false;
    
    // grab michel hits
    auto const& michel_hits = cluster._michel;

    // clear the map that links hit index to cluster they have been saved to
    _hit_cluster_map.clear();

    // clear this michel's cluster of photons
    cluster._michel._photon_clus_v.clear();
    cluster._michel._electron_hit_idx_v.clear();

    // keep track of which hit is closest to the Michel start point
    // this will be our electron cluster
    auto const& boundaryHit = cluster._hits[  cluster._boundary ];
    size_t electron_idx = 0;
    double dmin = 1030*1030;

    // loop trhough all hits in michel
    for (size_t i=0; i < michel_hits.size(); i++){

      auto const& pt1 = michel_hits[i];

      // search for hit closest to Michel start point
      if (pt1.SqDist( boundaryHit ) < dmin){
	electron_idx = i;
	dmin = pt1.SqDist( boundaryHit );
      }

      for (size_t j=i+1; j < michel_hits.size(); j++){

	auto const& pt2 = michel_hits[j];

	//std::cout << "exploring pair: " << i << " and " << j << std::endl;

	// if close enough
	if (pt1.SqDist(pt2) < _d_max) {

	  //std::cout << "they are close! merge.." << std::endl;

	  // if both hits already recorded -> in separate clusters
	  // merge their clusters
	  if ( ( _hit_cluster_map.find( i ) != _hit_cluster_map.end() ) and
	       ( _hit_cluster_map.find( j ) != _hit_cluster_map.end() ) ) {

	    //std::cout << "found hits " << i << " and " << j << " to belong to same cluster and already have different assigned clusters.." << std::endl;

	    // make sure these are different clusters
	    if ( _hit_cluster_map[i] != _hit_cluster_map[j] ) {
	      
	      // which clusters need to be merged?
	      auto idx1 = _hit_cluster_map[i];
	      auto idx2 = _hit_cluster_map[j];

	      //std::cout << "hits associated w/ clusters " << idx1 << " and " << idx2 << std::endl;

	      auto cl1 = cluster._michel._photon_clus_v[ idx1 ];
	      auto cl2 = cluster._michel._photon_clus_v[ idx2 ];

	      //std::cout << "w/ size " << cl1.size() << " and " << cl2.size() << std::endl;

	      // create a merged cluster
	      std::vector<size_t> clnew;
	      for (auto const& idx : cl1)
		clnew.push_back( idx );
	      for (auto const& idx : cl2)
		clnew.push_back( idx );

	      // for all indices in these clusters, replace the map that points
	      // to their cluster position in _michel._photon_clus_v
	      // with what will become the next cluster
	      for (auto const& idx : cl1)
		_hit_cluster_map[idx] = cluster._michel._photon_clus_v.size();
	      for (auto const& idx : cl2)
		_hit_cluster_map[idx] = cluster._michel._photon_clus_v.size();

	      //std::cout << "remove indices " << idx1 << " and " << idx2 << " from vec of size " << cluster._michel._photon_clus_v.size() << std::endl;

	      // erase the two clusters and replace with the merged version
	      cluster._michel._photon_clus_v[ idx1 ] = {};
	      cluster._michel._photon_clus_v[ idx2 ] = {};
	      // start with the largest of the two
	      /*
	      if (idx1 > idx2){
		cluster._michel._photon_clus_v.erase( cluster._michel._photon_clus_v.begin() + idx1 );
		cluster._michel._photon_clus_v.erase( cluster._michel._photon_clus_v.begin() + idx2 );
	      }
	      else{
		cluster._michel._photon_clus_v.erase( cluster._michel._photon_clus_v.begin() + idx2 );
		cluster._michel._photon_clus_v.erase( cluster._michel._photon_clus_v.begin() + idx1 );
	      }
	      */
	      // add new cluster
	      cluster._michel._photon_clus_v.push_back( clnew );
	      // for all indices in clnew -> add them to map
	      for (auto const& idx : clnew)
		_hit_cluster_map[idx] = cluster._michel._photon_clus_v.size() - 1;
	      
	    }// if both hits belong to different clusters
	    else{
	      if(_verbosity >= msg::kWARNING) {
		std::stringstream ss;
		ss << "\n\t\t Hits belong to same cluster! Error " << std::endl;
		Print(msg::kERROR,__FUNCTION__,ss.str());
	      }
	    }// if both belong to same cluster -> some thing went wrong
	  }// if both hits already assigned to a cluster

	  // is hit i recorded in the map?
	  else if ( _hit_cluster_map.find( i ) != _hit_cluster_map.end() ){
	    //std::cout << " hit i already in a cluster" << std::endl;
	    //std::cout << " hit i goes to cluster " << _hit_cluster_map[i] << std::endl;
	    //std::cout << " total # of clusters : " << cluster._michel._photon_clus_v.size() << std::endl;
	    cluster._michel._photon_clus_v[ _hit_cluster_map[ i ] ].push_back( j );
	    _hit_cluster_map[j] = _hit_cluster_map[i];
	  }
	  // is hit j recorded in the map?
	  else if ( _hit_cluster_map.find( j ) != _hit_cluster_map.end() ){
	    //std::cout << " hit j already in a cluster" << std::endl;
	      cluster._michel._photon_clus_v[ _hit_cluster_map[ j ] ].push_back( i );
	      _hit_cluster_map[i] = _hit_cluster_map[j];
	  }
	  else{
	    //std::cout << " neither hit in a cluster" << std::endl;
	    std::vector<size_t> new_clus = {i,j};
	    cluster._michel._photon_clus_v.push_back( new_clus );
	    _hit_cluster_map[i] = cluster._michel._photon_clus_v.size() - 1;
	    _hit_cluster_map[j] = cluster._michel._photon_clus_v.size() - 1;
	  }
	  
	}// if points are close
	
      }// 2nd loop throguh hits
    }// 1st loop through hits


    // are there isolated hits? if so add them as their own cluster
    for (size_t i=0; i < michel_hits.size(); i++){

      // if hit in no cluster
      if ( _hit_cluster_map.find( i ) == _hit_cluster_map.end() ){
	// add to a new 1-hit cluster
	std::vector<size_t> newclus = {i};
	cluster._michel._photon_clus_v.push_back( newclus );
	_hit_cluster_map[i] = cluster._michel._photon_clus_v.size() - 1 ;
	
      }// if hit in no cluster
    }// for all hits

    // find the cluster associated with the electron index
    auto clus_idx = _hit_cluster_map[electron_idx];

    //if there are only 1 or 2 hits in the Michel electron portion -> remove
    if ( cluster._michel._photon_clus_v.at( clus_idx ).size() <= 2)
      return false;

    // save electron cluster
    cluster._michel._electron_hit_idx_v = cluster._michel._photon_clus_v.at( clus_idx );
    cluster._michel._photon_clus_v.erase( cluster._michel._photon_clus_v.begin() + clus_idx );

    // erase any photon clusters of 0 size
    auto photon_clusters = cluster._michel._photon_clus_v;
    cluster._michel._photon_clus_v.clear();
    for (auto const& photon : photon_clusters){
      if (photon.size() != 0)
	cluster._michel._photon_clus_v.push_back( photon );
    }

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
