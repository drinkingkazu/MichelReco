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

    std::cout << "remove bad photons " << std::endl;
    
    ClusterVectorCalculator _clusterCalc;

    // loop through all photon clusters reconstructed
    // if cluster has too many hits, too much charge
    // or is too linear -> remove as probable
    // contamination from nearby muons

    auto const& michel = cluster._michel;

    // get all photon clusters
    auto const& photon_clus_v = michel._photon_clus_v;

    std::cout << "There are " << photon_clus_v.size() << " photon clusters..." << std::endl;

    // loop through all photon clusters and save only ones that
    // pass quality cuts
    std::vector< std::vector<size_t> > good_photon_clus_v;
    
    for (size_t i=0; i < photon_clus_v.size(); i++){

      auto const& photon_hit_v = photon_clus_v.at(i);

      // if too many hits -> remove
      if (photon_hit_v.size() > 20){
	std::cout << "removing photon cluster w/ " << photon_hit_v.size() << " hits" << std::endl;
	continue;
      }

      // if too many hits -> remove
      if (photon_hit_v.size() < 2){
	std::cout << "removing photon cluster w/ " << photon_hit_v.size() << " hits" << std::endl;
	continue;
      }	
      
      // very few hits? save
      if (photon_hit_v.size() < 5){
	good_photon_clus_v.push_back( photon_hit_v );
	std::cout << "less than 5 hits...keep this event." << std::endl;
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

      std::cout << "measured a linearity of " << linearity << std::endl;
      
      // too high? continue
      if (linearity < 0)
	linearity *= -1;
      if ( linearity > _max_linearity){
	std::cout << "remove photon cluster due to high linearity : " << linearity << std::endl;
	continue;
      }
      
      // made it this far -> add to good clusters
      good_photon_clus_v.push_back( photon_hit_v );
      
    }// for this photon
      
    cluster._michel._photon_clus_v = good_photon_clus_v;

    return true;
  }

}

#endif
