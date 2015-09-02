#ifndef MICHELRECO_CXX
#define MICHELRECO_CXX

#include "MichelReco.h"
#include "MichelException.h"
namespace michel {

  //-----------------------------------------------------------------
  MichelReco::MichelReco()
  //-----------------------------------------------------------------
    : _d_cutoff           ( 6.0 ) //Used to be 3.6
    , _min_nhits          ( 4   ) //Used to be 25
    , _verbosity          ( msg::kNORMAL )

    , _alg_merge          ( nullptr )
      
    , _alg_v        ( kAlgoTypeMax, nullptr )
    , _alg_time_v   ( kAlgoTypeMax, 0.      )
    , _alg_ctr_v    ( kAlgoTypeMax, 0       )
  {}
  
  //-----------------------------------------------------------------
  void MichelReco::Append(const std::vector<michel::HitPt>& hit_v)
  //-----------------------------------------------------------------
  {
    if(hit_v.size() < _min_nhits) return;
    MichelCluster cluster(_min_nhits, _d_cutoff);
    cluster.SetVerbosity(_verbosity);
    cluster.SetHits(hit_v);
    if(cluster._hits.size() < _min_nhits) return;
    _input_v.emplace_back(cluster);
  }

  //---------------------------------------------------------
  void MichelReco::Append(std::vector<michel::HitPt>&& hit_v)
  //---------------------------------------------------------
  {
    if(hit_v.size() < _min_nhits) return;
    MichelCluster cluster(std::move(hit_v), _min_nhits, _d_cutoff);
    cluster.SetVerbosity(_verbosity);
    if(cluster._hits.size() < _min_nhits) return;
    _input_v.emplace_back(cluster);
  }

  //---------------------------------------------------------------------------
  void MichelReco::RegisterAllHits(const std::vector<michel::HitPt>& all_hit_v)
  //---------------------------------------------------------------------------
  {
    if(all_hit_v.empty()){
      std::cerr << "\033[93m[ERROR]\033[00m cannot have hit & marker list w/ different length..."
		<< std::endl;
      throw MichelException();
    }
    _all_hit_v = all_hit_v;
  }

  //-----------------------------------------------------------------
  void MichelReco::RegisterAllHits(std::vector<michel::HitPt>&& all_hit_v)
  //-----------------------------------------------------------------
  {
    if(all_hit_v.empty()){
      std::cerr << "\033[93m[ERROR]\033[00m cannot have hit & marker list w/ different length..."
		<< std::endl;
      throw MichelException();
    }
    std::swap(_all_hit_v,all_hit_v);
  }
  
  //-----------------------------------------------------------------
  void MichelReco::AddAlgo(BaseMichelAlgo* algo)
  //-----------------------------------------------------------------
  {
    _alg_v.push_back(algo);
    _alg_time_v.push_back(0.);
    _alg_ctr_v.push_back(0);
  }
  
  //-----------------------------------------------------------------
  void MichelReco::AddAna(MichelAnaBase* ana)
  //-----------------------------------------------------------------
  {
    _ana_v.push_back(ana);
  }

  //-----------------------------------------------------------------
  void MichelReco::Initialize()
  //-----------------------------------------------------------------
  {
    for(auto& ana : _ana_v) {
      ana->SetVerbosity(_verbosity);
      ana->Initialize();
    }

    for(size_t i=0; i<_alg_v.size(); ++i) {

      if(_alg_v[i]) _alg_v[i]->SetVerbosity(_verbosity);
      
      if(i == kClusterMerger) continue;

      if(_alg_v[i]) continue;

      std::cerr << "\033[93m[ERROR]\033[00m "
		<< "Algorithm type : " << i
		<< " not provided! " << std::endl;
      throw MichelException();
    }
  }
  
  //-----------------------------------------------------------------
  void MichelReco::Process()
  //-----------------------------------------------------------------
  {
   	
    if(_verbosity <= msg::kDEBUG)
      std::cout << "<<Process>> start!" << std::endl;

    // If nothing to be done, return
    if(_input_v.empty()) return;

    // make sure hit vectors provided & cluster list (another hit arrays) are consistent
    _used_hit_marker_v.resize(_all_hit_v.size(),false);
    for(size_t i=0; i<_used_hit_marker_v.size(); ++i) {
      if(_all_hit_v[i]._pl != 2) //check plane
	_used_hit_marker_v[i] = true;
      else
	_used_hit_marker_v[i] = false;
    }

    if(_verbosity <= msg::kDEBUG)

      std::cout << "<<Process>> running kClusterMerge..." << std::endl;
    
    //
    // Step 1 ... merge clusters
    //
    if(!_alg_merge)

      _output_v = _input_v;

    else {
      _watch.Start();
      _output_v = _alg_merge->Merge(_input_v);
      _alg_time_v [kClusterMerger] += _watch.RealTime();
      _alg_ctr_v  [kClusterMerger] += _input_v.size();
    }
    
    // Update "used hits" list
    for(auto const& cluster : _output_v) {
     
      for(auto const& hit_pt : cluster._hits) {

	if(hit_pt._id >= _used_hit_marker_v.size()) {

	  std::cout << "\033[93m[ERROR]\033[00m "
		    << "Found a hit index out of range! "
		    << hit_pt._id << " out of " << _used_hit_marker_v.size()
		    << std::endl;
	  throw MichelException();
	}
	_used_hit_marker_v[hit_pt._id] = true;
      }
    }
    
    if(_verbosity <= msg::kDEBUG)

      std::cout << "<<Process>> running kBoundayFinder..." << std::endl;

    for(auto& cluster : _output_v ){
      // start going through algorithms and executing
      // them consecutively
      for (size_t n=0; n < _alg_v.size(); n++){
	if(_verbosity <= msg::kDEBUG) { std::cout << "<<Process>> running algo : " << _alg_v[n]->Name() << std::endl; }
	_watch.Start();
	_alg_v[n]->ProcessCluster(cluster,_all_hit_v);
	_alg_time_v[n] += _watch.RealTime();
	_alg_ctr_v[n] += 1;
	// if the cluster is now a default instance of MichelCluster
	// then don't bother to move on to the next algorithms
	if (cluster._has_michel == false)
	  break;
      }// looping through algorithms
    }// for all clusters


    //
    // Finally call analyze
    //
    for(auto& ana : _ana_v){
      // firs set the event id for the ana
      ana->SetEventID(_id);
      ana->Analyze(_input_v,_output_v);
    }

    if(_verbosity <= msg::kDEBUG)

      std::cout << "<<Process>> end!" << std::endl;
  }
  
  //-----------------------------------------------------------------
  void MichelReco::EventReset()
  //-----------------------------------------------------------------
  {
    for(auto& alg : _alg_v) if(alg) alg->EventReset();
    for(auto& ana : _ana_v) ana->EventReset();
    _input_v.clear();
    _output_v.clear();
  }

  //-----------------------------------------------------------------
  void MichelReco::Finalize(TFile *fout)
  //-----------------------------------------------------------------
  {

    for(auto& ana : _ana_v) ana->Finalize(fout);

    // loop through algos and evaluate time-performance
    std::cout << std::endl
	      << "=================== Time Report =====================" << std::endl;
    for(size_t n=0; n < _alg_v.size(); n++){
      double alg_time = _alg_time_v[n]/((double)_alg_ctr_v[n]);
      std::cout <<  _alg_v[n]->Name() << " Algo Time: " << alg_time*1.e6     << " [us/cluster]" << std::endl;
    }
      
    std::cout << "=====================================================" << std::endl
	      << std::endl;
  }
  

}
#endif
