#ifndef COSMICANA_CXX
#define COSMICANA_CXX

#include "CosmicAna.h"

namespace michel {

  /// Initialize
  void CosmicAna::Initialize()
  {
    std::cout << "Initializing.. tree...\n";
    _out_tree = new TTree("out_tree","aho_tree");

    _out_tree->Branch("_michel_clustered_charge", "std::vector<double>", &_michel_clustered_charge);
    _out_tree->Branch("_michel_n_hits"          , "std::vector<double>", &_michel_n_hits          );
    _out_tree->Branch("_number_of_clusters"     , &_number_of_clusters , "_number_of_clusters/I");

    _out_tree->Branch("_boundary", "std::vector<int>", &_boundary);

    _out_tree->Branch("_chi_at_boundary", "std::vector<double>", &_chi_at_boundary);
    _out_tree->Branch("_mean_chi"       , "std::vector<double>", &_mean_chi);
    _out_tree->Branch("_rms_chi"        , "std::vector<double>", &_rms_chi);
    _out_tree->Branch("_lowest_chi"     , "std::vector<double>", &_lowest_chi);
    
    _out_tree->Branch("_Z", "std::vector<double>" , &_Z);
    _out_tree->Branch("_X", "std::vector<double>" , &_X);
    
    _out_tree->Branch("_michel_Z", "std::vector<double>" , &_michel_Z);
    _out_tree->Branch("_michel_X", "std::vector<double>" , &_michel_X);

    _out_tree->Branch("_has_michel", "std::vector<bool>" , &_has_michel);
    
    
  }
  
  /// Analyze
  void CosmicAna::Analyze(const MichelClusterArray& input_cluster_v,
			  const MichelClusterArray& output_cluster_v)
  {
    clear_all();
    std::vector<double> EMPTYVEC = {-500};
    
    // Loop over all output clusters, put X,Z into _Z,_X
    _Z.reserve( output_cluster_v.size() * 4 * 25);
    _X.reserve( output_cluster_v.size() * 4 * 25);

    _michel_clustered_charge.reserve(output_cluster_v.size());
    _michel_n_hits.reserve(output_cluster_v.size());

    _boundary.reserve(output_cluster_v.size());

    _chi_at_boundary.reserve(output_cluster_v.size());
    _mean_chi.reserve(output_cluster_v.size());
    _rms_chi.reserve(output_cluster_v.size());
    _lowest_chi.reserve(output_cluster_v.size());
    
    _michel_Z.reserve(output_cluster_v.size() * 4 * 25);
    _michel_X.reserve(output_cluster_v.size() * 4 * 25);
    
    _has_michel.reserve(output_cluster_v.size());
    
    
    std::vector<double> covariance_in_largest_cluster;
    covariance_in_largest_cluster.reserve( output_cluster_v.size() * 4 * 25);
    
    //get the number of clusters
    int number_of_good_clusters = 0;
    
    //fill cluster-wise info
    for(const auto& out : output_cluster_v) {
      
      if(out._chi2_v.size() == 0)
	continue;

      number_of_good_clusters++;      

      covariance_in_largest_cluster.clear();

      for(const auto& c : out._chi2_v)
	covariance_in_largest_cluster.push_back(std::abs(c));
      
      //parse the hits...
      for(const auto& hit : out._hits)  {
	_Z.push_back(hit._w); _X.push_back(hit._t);
      }

      //get the boundary
      auto boundary = (int)out._boundary;
      
      //parse the chi^2
      _boundary.push_back       ( boundary                                   );
      _chi_at_boundary.push_back( covariance_in_largest_cluster.at(boundary) );
      _mean_chi.push_back       ( get_mean   (covariance_in_largest_cluster) );
      _rms_chi.push_back        ( get_rms    (covariance_in_largest_cluster) );
      _lowest_chi.push_back     ( get_lowest (covariance_in_largest_cluster) );
      
      
      //if there is a michel...
      if(out._michel.size()) {
	_has_michel.push_back(true);
	auto total_charge = double{0.0};
	for(const auto& mhit : out._michel) {
	  _michel_Z.push_back(mhit._w); 
	  _michel_X.push_back(mhit._t);
	  total_charge += mhit._q;
	} 
	_michel_clustered_charge.push_back(total_charge);
	_michel_n_hits.push_back(out._michel.size());	
      }
      else {
	_michel_clustered_charge = EMPTYVEC;
	_michel_n_hits           = EMPTYVEC;
	_michel_Z                = EMPTYVEC;
	_michel_X                = EMPTYVEC;
	_has_michel.push_back(false);
      }
      

		
      
    }
    _number_of_clusters = number_of_good_clusters;
    
    _out_tree->Fill();
    

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
    
    _michel_clustered_charge.clear();
    _michel_n_hits.clear();
    _number_of_clusters = 0;

    _boundary.clear();

    _chi_at_boundary.clear();
    _mean_chi.clear();
    _rms_chi.clear();
    _lowest_chi.clear();
    
    _Z.clear();
    _X.clear();
    
    _michel_Z.clear();
    _michel_X.clear();
    
    
  }
  
}
#endif

