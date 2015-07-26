#ifndef AHOANA_CXX
#define AHOANA_CXX

#include "AhoAna.h"

namespace michel {

  /// Initialize
  void AhoAna::Initialize()
  {
    std::cout << "Initializing.. tree...\n";
    _out_tree = new TTree("out_tree","aho_tree");
    _out_tree->Branch("_largest_cluster_charge",&_largest_cluster_charge,"_largest_cluster_charge/D");
    _out_tree->Branch("_n_hits_in_largest_cluster",&_n_hits_in_largest_cluster,"_n_hits in_largest_cluster/D");
    _out_tree->Branch("_n_hits_in_largest_cluster_michel",&_n_hits_in_largest_cluster_michel,"_n_hits_in_largest_cluster_michel/D");
    
  }
  
  /// Analyze
  void AhoAna::Analyze(const MichelClusterArray& input_cluster_v,
		       const MichelClusterArray& output_cluster_v)
  {
    bool lets_get_physical = false;
    int size = 0;
    
    for(auto& output : output_cluster_v) {
      if(output._hits.size() > size) {
	if(output._michel.size()) {
	  _largest_cluster_charge           = output._michel._charge;
	  _n_hits_in_largest_cluster        = output._hits.size();
	  _n_hits_in_largest_cluster_michel = output._michel.size();
	  size = output._hits.size();
	  lets_get_physical = true;
	}
      }
    }
    
    //No real check on anything just fill it whatever
    if(lets_get_physical)
      _out_tree->Fill();
    
  }
  
  /// Event Reset
  void AhoAna::EventReset()
  {
    
  }
  
  /// Finalize
  void AhoAna::Finalize(TFile* fout)
  {
    _out_tree->Write();
    //fout->Write();
  }
  
}
#endif

