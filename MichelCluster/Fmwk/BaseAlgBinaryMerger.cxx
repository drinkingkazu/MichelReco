#ifndef BASEALGBINARYMERGER_CXX
#define BASEALGBINARYMERGER_CXX

#include "BaseAlgBinaryMerger.h"
#include "CMergeBookKeeper.h"
#include <map>

namespace michel {

  MichelClusterArray BaseAlgBinaryMerger::Merge(const MichelClusterArray& input_v)
  {

    MichelClusterArray result_v = input_v;

    if(input_v.size()<2) return result_v;
    
    while(_recursive) {
      // Instantiate a book keeper
      CMergeBookKeeper bk(result_v.size());

      // Compute the priority order & store in an ordered map
      std::multimap<double,size_t> score_index_m;
      for(size_t cluster_index=0; cluster_index<result_v.size(); ++cluster_index) {

	auto const& cluster = result_v[cluster_index];

	auto score = Priority(cluster);

	if( score <= 0 ) continue;

	score_index_m.emplace( 1./Priority(cluster), cluster_index );
      }

      // Loop over possible combinations and try merging
      for(size_t a_index=0; a_index<score_index_m.size(); ++a_index) {

	// Find cluster a's index
	auto score_iter_a = score_index_m.begin();
	std::advance(score_iter_a,a_index);
	// Obtain cluster a's reference
	auto const& cluster_a = result_v[(*score_iter_a).second];

	for(size_t b_index=(a_index+1); b_index<score_index_m.size(); ++b_index) {

	  // Find cluster b's index
	  auto score_iter_b = score_index_m.begin();
	  std::advance(score_iter_b,b_index);
	  // Obtain cluster b's reference
	  auto const& cluster_b = result_v[(*score_iter_b).second];

	  // Inspect to be merged or not
	  if( Merge(cluster_a,cluster_b) )

	    bk.Merge( (*score_iter_a).second, (*score_iter_b).second );
		      
	}	  
      }
      
      // If nothing merged, break
      auto const& bk_result = bk.GetResult();

      if(bk_result.size() == result_v.size()) break;
      else {
	// Merge ... prepare a temporary vector container
	MichelClusterArray tmp_result_v;
	tmp_result_v.reserve(bk_result.size());

	// Loop over merged cluster sets
	for(auto const& index_set : bk_result) {

	  // Construct a hit vector to be made for a cluster
	  std::vector<HitPt> hits;
	  double d_cutoff=0;
	  size_t min_nhits=0;
	  // prepare list of cluster indices that are being merged together
	  std::vector<unsigned short> input_clus_idx_v;
	  // Loop over index numbers of associated hits
	  for(auto const& index : index_set) {
	    auto const& cluster = result_v[index];
	    if(!min_nhits) {
	      min_nhits = cluster._min_nhits;
	      d_cutoff  = cluster._d_cutoff;
	    }
	    hits.reserve(hits.size() + cluster._hits.size());
	    for(auto const& h : cluster._hits)
	      hits.push_back(h);
	    // get input cluster index list from this MichelCluster
	    auto clus_idx_v = result_v[index].getInputClusterIndex_v();
	    for (auto& clus : clus_idx_v)
	      input_clus_idx_v.push_back(clus);
	  }
	  //MichelCluster merged(std::move(hits), min_nhits,d_cutoff);
	  if (hits.size() >= 3){
	    MichelCluster merged(std::move(hits), 3, d_cutoff);
	    // save the index set as this MichelCluster's list of input clusters
	    merged.setInputClusterIndex_v(input_clus_idx_v);
	    tmp_result_v.emplace_back(merged);
	  }
	}
	result_v = tmp_result_v;
      }
    }



    return result_v;
  }

}
#endif
