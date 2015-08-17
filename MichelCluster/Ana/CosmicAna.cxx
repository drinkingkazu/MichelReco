#ifndef COSMICANA_CXX
#define COSMICANA_CXX

#include "CosmicAna.h"

namespace michel {

  /// Initialize
  void CosmicAna::Initialize()
  {
    std::cout << "Initializing.. tree...\n";
    _out_tree = new TTree("out_tree","aho_tree");

    _out_tree->Branch("_michel_clustered_charge", &_michel_clustered_charge, "_michel_clustered_charge/D");
    _out_tree->Branch("_michel_n_hits"          , &_michel_n_hits, "_michel_n_hits/D"          );
    //_out_tree->Branch("_number_of_clusters"     , &_number_of_clusters , "_number_of_clusters/I");

    _out_tree->Branch("_boundary", &_boundary, "_boundary/I");

    _out_tree->Branch("_chi_at_boundary", &_chi_at_boundary, "_chi_at_boundary/D");
    _out_tree->Branch("_mean_chi"       , &_mean_chi       , "_mean_chi/D");
    _out_tree->Branch("_rms_chi"        , &_rms_chi        , "_rms_chi/D");
    _out_tree->Branch("_lowest_chi"     , &_lowest_chi     , "_lowest_chi/D");
    
    _out_tree->Branch("_Z", "std::vector<double>" , &_Z);
    _out_tree->Branch("_X", "std::vector<double>" , &_X);
    
    _out_tree->Branch("_michel_Z", "std::vector<double>" , &_michel_Z);
    _out_tree->Branch("_michel_X", "std::vector<double>" , &_michel_X);

    _out_tree->Branch("_q_v",    "std::vector<double>" , &_q_v);

    _out_tree->Branch("_t_q_v",    "std::vector<double>" , &_t_q_v);
    _out_tree->Branch("_t_dqds_v", "std::vector<double>" , &_t_dqds_v);
    _out_tree->Branch("_chi_v",    "std::vector<double>" , &_chi_v);

    _out_tree->Branch("_s_v",    "std::vector<double>" , &_s_v);

    _out_tree->Branch("_has_michel", &_has_michel, "_has_michel/O");
    
    
  }
  
  /// Analyze
  void CosmicAna::Analyze(const MichelClusterArray& input_cluster_v,
			  const MichelClusterArray& output_cluster_v)
  {
    clear_all();
    std::vector<double> EMPTYVEC = {};
    
    // Loop over all output clusters, put X,Z into _Z,_X
    _Z.reserve( output_cluster_v.size() * 4 * 25);
    _X.reserve( output_cluster_v.size() * 4 * 25);

    _michel_Z.reserve(output_cluster_v.size() * 4 * 25);
    _michel_X.reserve(output_cluster_v.size() * 4 * 25);
    _t_q_v.reserve   (1000);
    _t_dqds_v.reserve(1000);
    _chi_v.reserve   (1000);
    
    std::vector<double> covariance_in_largest_cluster;
    covariance_in_largest_cluster.reserve( output_cluster_v.size() * 4 * 25);
    
    //get the number of clusters
    int number_of_good_clusters = 0;
    
    //fill cluster-wise info
    for(const auto& out : output_cluster_v) {

      clear_all(); 
      covariance_in_largest_cluster.clear();

      if(out._chi2_v.size() == 0)
	continue;
      
      for(const auto& c : out._chi2_v) {
	covariance_in_largest_cluster.push_back( abs(c) );
	_chi_v.push_back                       ( abs(c) );
      }
      
      //parse the hits...
      for(const auto& hit : out._hits)  {
	_Z.push_back(hit._w); _X.push_back(hit._t); _q_v.push_back(hit._q);
      }
      
      _t_q_v    = out._t_mean_v;
      _t_dqds_v = out._t_dqds_v;
      _s_v      = out._s_v;

      //get the boundary
      auto boundary = (int)out._boundary;
      
      //parse the chi^2
      _boundary         =  boundary                                   ;
      _chi_at_boundary  =  covariance_in_largest_cluster.at(boundary) ;
      _mean_chi         =  get_mean   (covariance_in_largest_cluster) ;
      _rms_chi          =  get_rms    (covariance_in_largest_cluster) ;
      _lowest_chi       =  get_lowest (covariance_in_largest_cluster) ;
      
      
      //if there is a michel...
      if(out._michel.size()) {
	
	_has_michel = true;
	auto total_charge = double{0.0};
	
	for(const auto& mhit : out._michel) {
	  _michel_Z.push_back(mhit._w); 
	  _michel_X.push_back(mhit._t);
	  total_charge += mhit._q;
	} 
	
	_michel_clustered_charge = total_charge;
	_michel_n_hits           = out._michel.size();

      }
      else {
	_michel_clustered_charge = -1;
	_michel_n_hits           = -1;
	_michel_Z                = EMPTYVEC;
	_michel_X                = EMPTYVEC;
	_has_michel              = false;
      }
		
      _out_tree->Fill();      
    }
    
  }
  
  /// Event Reset
  void CosmicAna::EventReset()
  {
    
  }
  
  /// Finalize
  void CosmicAna::Finalize(TFile* fout)
  {
    _out_tree->Write();
    //fout->Write();
  }



  double CosmicAna::get_mean(const std::vector<double>& data) {

    if(data.size() == 0){ std::cout << "You have me nill to mean\n"; throw MichelException(); }
    
    double result = 0.0;

    for(const auto& d : data) 
      result += d;
    
    
    return (result / ((double)data.size()));

    
  }
  double CosmicAna::get_rms(const std::vector<double>& data){

    if(data.size() == 0){ std::cout << "You have me nill to stdev\n"; throw MichelException(); }

    double result = 0.0;
    auto    avg   = get_mean(data);
    for(const auto& d: data)
      result += (d - avg)*(d - avg);
    
    return sqrt(result/((double)data.size()));

    
  }
  double CosmicAna::get_lowest(const std::vector<double>& data){
    
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
  
  void CosmicAna::clear_all() {
    
    _michel_clustered_charge = -1;
    _michel_n_hits           = -1;
    _number_of_clusters      = -1;

    _boundary = -1;

    _chi_at_boundary = -1;
    _mean_chi        = -1;
    _rms_chi         = -1;
    _lowest_chi      = -1;
    
    _Z.clear();
    _X.clear();
    
    _michel_Z.clear();
    _michel_X.clear();
    
    _q_v.clear();
    _t_q_v.clear();

    _t_dqds_v.clear();
    _chi_v.clear();
    
    _s_v.clear();
  }
  
}
#endif

