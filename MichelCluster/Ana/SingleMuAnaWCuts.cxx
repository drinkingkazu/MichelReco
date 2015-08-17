#ifndef SINGLEMUANAWCUTS_CXX
#define SINGLEMUANAWCUTS_CXX

#include "SingleMuAnaWCuts.h"

namespace michel {

  /// Initialize
  void SingleMuAnaWCuts::Initialize()
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
    _out_tree->Branch("mean_chi",   &mean_chi, "mean_chi/D");
    _out_tree->Branch("rms_chi",    &rms_chi, "rms_chi/D");
    _out_tree->Branch("lowest_chi", &lowest_chi, "lowest_chi/D");
	
    
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
  void SingleMuAnaWCuts::Analyze(const MichelClusterArray& input_cluster_v,
			    const MichelClusterArray& output_cluster_v)
  {
    
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


    chi_at_boundary = -1.0;
    mean_chi        = -1.0;
    rms_chi         = -1.0;
    lowest_chi      = -1.0;



    
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
	
	for(const auto& c : output._chi2_v)
	  covariance_in_largest_cluster.push_back(abs(c));

	
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

    if(covariance_in_largest_cluster.size() > 0) {
      chi_at_boundary = covariance_in_largest_cluster.at(boundary);
      mean_chi        = get_mean(covariance_in_largest_cluster);
      rms_chi         = get_rms(covariance_in_largest_cluster);
      lowest_chi      = get_lowest(covariance_in_largest_cluster);
    }


    //From ipython notebooks:
    // S = SIGNAL.query('_largest_cluster_charge > 0')
    //   S = S.query('_number_of_clusters >= 1')
    //   S = S.query('_n_hits_in_largest_cluster_michel >= 10')
    //   S = S.query('lowest_chi < 0.28')
    //   S = S.query('chi_at_boundary < 0.68')
    //   S = S.query('mean_chi > 0.98')

    if(_largest_cluster_charge > 0 && 
       _number_of_clusters >= 1 &&
       _n_hits_in_largest_cluster_michel >= 10 &&
       lowest_chi < 0.28 &&
       chi_at_boundary < 0.68 &&
       mean_chi > 0.98 &&
       _largest_cluster_charge < 8929)
      
      _out_tree->Fill();
    
    
    
  }
  
  /// Event Reset
  void SingleMuAnaWCuts::EventReset()
  {
    
  }
  
  /// Finalize
  void SingleMuAnaWCuts::Finalize(TFile* fout)
  {
    _out_tree->Write();
    //fout->Write();
  }



  double SingleMuAnaWCuts::get_mean(const std::vector<double>& data) {

    if(data.size() == 0){ std::cout << "You have me nill to mean\n"; throw MichelException(); }
    
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
    
    
    return (result / ((double)data.size()));

    
  }
  double SingleMuAnaWCuts::get_rms(const std::vector<double>& data){

    if(data.size() == 0){ std::cout << "You have me nill to stdev\n"; throw MichelException(); }

    double result = 0.0;
    auto    avg   = get_mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));

    
  }
  double SingleMuAnaWCuts::get_lowest(const std::vector<double>& data){
    
    //get lowest dqds
    auto the_min = double{999999.0};
    size_t idx= 1;
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] < the_min) {
	the_min = data[i]; idx = i;
      }
    }
    
    return the_min;

    
  }
		    
  
}
#endif

