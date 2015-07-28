#ifndef AHOANA_CXX
#define AHOANA_CXX

#include "AhoAna.h"

namespace michel {

  /// Initialize
  void AhoAna::Initialize()
  {
    std::cout << "Initializing.. tree...\n";
    _out_tree = new TTree("out_tree","aho_tree");
    _out_tree->Branch("_largest_cluster_charge",          &_largest_cluster_charge,          "_largest_cluster_charge/D");
    _out_tree->Branch("_n_hits_in_largest_cluster",       &_n_hits_in_largest_cluster,       "_n_hits in_largest_cluster/D");
    _out_tree->Branch("_n_hits_in_largest_cluster_michel",&_n_hits_in_largest_cluster_michel,"_n_hits_in_largest_cluster_michel/D");

    _out_tree->Branch("_number_of_clusters",&_number_of_clusters,"_number_of_clusters/I");

    _out_tree->Branch("boundary", &boundary, "boundary/I");
    _out_tree->Branch("IMAX", &IMAX, "IMAX/I");
    _out_tree->Branch("IMIN", &IMIN, "IMIN/I");

    _out_tree->Branch("chi_at_boundary", &chi_at_boundary, "chi_at_boundary/D");
    
    _out_tree->Branch("_Z", "std::vector<double>" , &_Z);
    _out_tree->Branch("_X", "std::vector<double>" , &_X);

    _out_tree->Branch("michel_Z", "std::vector<double>" , &michel_Z);
    _out_tree->Branch("michel_X", "std::vector<double>" , &michel_X);
    
    _out_tree->Branch("charge_in_largest_cluster",          "std::vector<double>", &charge_in_largest_cluster);
    _out_tree->Branch("truncated_charge_in_largest_cluster","std::vector<double>", &truncated_charge_in_largest_cluster);
    _out_tree->Branch("truncated_dqds_in_largest_cluster",  "std::vector<double>", &truncated_dqds_in_largest_cluster);

    _out_tree->Branch("covariance_in_largest_cluster",  "std::vector<double>", &covariance_in_largest_cluster);
    
    _out_tree->Branch("s","std::vector<double>", &s);
    
  }
  
  /// Analyze
  void AhoAna::Analyze(const MichelClusterArray& input_cluster_v,
		       const MichelClusterArray& output_cluster_v)
  {
    
    //Write out this event to the TTree no matter what, whether I see the michel or not
    
    for(const auto& out : output_cluster_v) {
      for(const auto& hit : out._hits)  {
	_Z.push_back(hit._w); _X.push_back(hit._t);
      }

      for(const auto& mhit : out._michel) {
	michel_Z.push_back(mhit._w); michel_X.push_back(mhit._t);
      }
    }
    
    
    bool lets_get_physical = false;
    int sizE = -1;
    _number_of_clusters = 0;
    
    for(auto& output : output_cluster_v) {
      _number_of_clusters++;  

      if((int)output._hits.size() > sizE) {
	charge_in_largest_cluster.clear();
	truncated_charge_in_largest_cluster.clear();
	truncated_dqds_in_largest_cluster.clear();
	covariance_in_largest_cluster.clear();
	s.clear();
	
	for(const auto& hits : output._hits)
	  charge_in_largest_cluster.push_back(hits._q);
	
	for(const auto& c : output._t_mean_v)
	  truncated_charge_in_largest_cluster.push_back(c);

	for(const auto& c : output._t_dqds_v)
	  truncated_dqds_in_largest_cluster.push_back(c);

	for(const auto& c : output._s_v)
	  s.push_back(c);
	
	std::cout << "{";
	for(const auto& c : output._chi2_v) {
	  covariance_in_largest_cluster.push_back(c);
	  std::cout << c << ",";
	  if(isnan(c)) { throw MichelException(); }
	}
	std::cout << "}";

	// truncated_charge_in_largest_cluster = output._t_mean_v;
	// truncated_dqds_in_largest_cluster   = output._t_dqds_v;
	// s                                   = output._s_v;
	sizE     = (int)output._hits.size();
	boundary = (int)output._boundary;

	if(output._michel.size()) {
	  _largest_cluster_charge           = output._michel._charge;
	  _n_hits_in_largest_cluster        = output._hits.size();
	  _n_hits_in_largest_cluster_michel = output._michel.size();
	  lets_get_physical = true;
	}
      }
    }
    
    //No real check on anything just fill it whatever

    _out_tree->Fill();
    
    //Clear variables
    _largest_cluster_charge           = -1.0;
    _n_hits_in_largest_cluster        = -1.0;
    _n_hits_in_largest_cluster_michel = -1.0;
    _Z.clear(); _X.clear();
    michel_Z.clear(); michel_X.clear();
    charge_in_largest_cluster.clear();
    truncated_charge_in_largest_cluster.clear();
    truncated_dqds_in_largest_cluster.clear();
    covariance_in_largest_cluster.clear();
    s.clear();
    boundary = -1;
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

