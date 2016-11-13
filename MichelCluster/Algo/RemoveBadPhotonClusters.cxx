#ifndef REMOVEBADPHOTONCLUSTERS_CXX
#define REMOVEBADPHOTONCLUSTERS_CXX

#include "RemoveBadPhotonClusters.h"
#include "Fmwk/ClusterVectorCalculator.h"

#include <sstream>

namespace michel{

  RemoveBadPhotonClusters::RemoveBadPhotonClusters()
  {
    _name    = "RemoveBadPhotonClusters";
    _max_linearity = 0.9;
  }

  bool RemoveBadPhotonClusters::ProcessCluster(MichelCluster& cluster,
					       const std::vector<HitPt>& hits)
  {

    if (hits.size() == 0) return false;
    
    ClusterVectorCalculator _clusterCalc;

    // loop through all photon clusters reconstructed
    // if cluster has too many hits, too much charge
    // or is too linear -> remove as probable
    // contamination from nearby muons

    auto const& michel = cluster._michel;

    // get all photon clusters
    auto const& photon_clus_v = michel._photon_clus_v;

    // if there are no phton clusters -> exit
    if (photon_clus_v.size() == 0)
      return true;

    // get michel start point location
    auto const& boundary = cluster._hits[ cluster._boundary ];

    // grab all hits and select those within 1.5 meters of the muon decay point
    // these will be used to check for conflicts
    std::vector<HitPt> hit_v;
    for (auto const& hit : hits){

      if (hit.SqDist(boundary) < 150*150)
	hit_v.push_back( hit );

    }// for all hits
      

    //std::cout << "There are " << photon_clus_v.size() << " photon clusters..." << std::endl;

    // loop through all photon clusters and save only ones that
    // pass quality cuts
    std::vector< std::vector<size_t> > good_photon_clus_v;
    
    for (size_t i=0; i < photon_clus_v.size(); i++){

      auto const& photon_hit_v = photon_clus_v.at(i);

      // if too many hits -> remove
      if (photon_hit_v.size() > 20){
	//std::cout << "removing photon cluster w/ " << photon_hit_v.size() << " hits" << std::endl;
	continue;
      }

      // if too many hits -> remove
      if (photon_hit_v.size() < 2){
	//std::cout << "removing photon cluster w/ " << photon_hit_v.size() << " hits" << std::endl;
	continue;
      }	
      
      // very few hits? save
      if (photon_hit_v.size() < 5){
	good_photon_clus_v.push_back( photon_hit_v );
	//std::cout << "less than 5 hits...keep this event." << std::endl;
	continue;
      }
      
      // measure linearity
      std::vector<double> photon_w_v;
      std::vector<double> photon_t_v;
      for (auto const& hit_idx : photon_hit_v){
	photon_w_v.push_back( michel.at( hit_idx )._w );
	photon_t_v.push_back( michel.at( hit_idx )._t );
      }
      double linearity = _clusterCalc.cov( photon_w_v, photon_t_v );

      //std::cout << "measured a linearity of " << linearity << std::endl;
      
      // too high? continue
      if (linearity < 0)
	linearity *= -1;
      if ( linearity > _max_linearity){
	//std::cout << "remove photon cluster due to high linearity : " << linearity << std::endl;
	continue;
      }

      // another check -> if far-away photons are surrounded by more charge
      // in small region surrounding those hits -> ignore because they
      // are probably from another muon
      // how many are nearby?
      int n_close = 0;
      // procedure :
      // 0) if photon cluster is close to Michel -> don't perform this check
      bool close = false;
      for (auto const& photon_hit_idx : photon_hit_v){
	auto const& photon_hit = michel.at( photon_hit_idx );
	if (photon_hit.Dist(boundary) < 20){
	  close = true;
	  break;
	}// if close
      }// for all photon hits
      if (close == false){
	// 1) find how many hits are within 1 cm of any hit in this cluster.
	//std::cout << "There are " << photon_hit_v.size() << " hits in this photon" << std::endl;
	for (auto const& hit : hit_v){
	  bool in_photon = false;
	  for (auto const& photon_hit_idx : photon_hit_v){
	    auto const& photon_hit = michel.at( photon_hit_idx );
	    // if hit index is already a photon hit, ignore
	    if (hit._id == photon_hit._id){
	      in_photon = true;
	      continue;
	    }
	  }// for all photon hits
	  // is this a photon hit? if so don't check
	  if (in_photon == true)
	    continue;
	  // else, loop through photon-htis to find distance
	  for (auto const& photon_hit_idx : photon_hit_v){
	    auto const& photon_hit = michel.at( photon_hit_idx );
	    // check distance
	    if (photon_hit.SqDist(hit) < 2.*2.){
	      n_close += 1;
	      // we already found that this hit is in the circle -> break loop over photon hits
	      break;
	    }// if hit is close
	  }// for all photon hits
	}// for all hits
	//std::cout << "there are " << n_close << " other hits narby" << std::endl;
      }//if cluster is not close
      
      // if there are hits nearby -> remove
      if (n_close > 0)
	continue;
      
      // made it this far -> add to good clusters
      //std::cout << "add this photon." << std::endl;
      good_photon_clus_v.push_back( photon_hit_v );

    }// for this photon
      
    cluster._michel._photon_clus_v = good_photon_clus_v;

    return true;
  }

}

#endif
