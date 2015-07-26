#ifndef MICHELRECO_CXX
#define MICHELRECO_CXX

#include "MichelReco.h"
#include "MichelException.h"
namespace michel {

  //-----------------------------------------------------------------
  MichelReco::MichelReco()
  //-----------------------------------------------------------------
    : _d_cutoff           ( 3.6 )
    , _min_nhits          ( 10  )
    , _verbosity          ( msg::kNORMAL )

    , _alg_merge          ( nullptr )
    , _alg_boundary       ( nullptr )
    , _alg_michel_id      ( nullptr )
    , _alg_michel_cluster ( nullptr )

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
    if(cluster._hits.empty()) return;
    _input_v.emplace_back(cluster);
  }

  //---------------------------------------------------------
  void MichelReco::Append(std::vector<michel::HitPt>&& hit_v)
  //---------------------------------------------------------
  {
    if(hit_v.size() < _min_nhits) return;
    MichelCluster cluster(std::move(hit_v), _min_nhits, _d_cutoff);
    cluster.SetVerbosity(_verbosity);
    if(cluster._hits.empty()) return;
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
  void MichelReco::SetAlgo(const AlgoType_t type,
			   BaseMichelAlgo* algo)
  //-----------------------------------------------------------------
  {
    switch(type) {
    case kClusterMerger:  _alg_merge          = (BaseAlgMerger*)algo;        break;
    case kBoundaryFinder: _alg_boundary       = (BaseAlgBoundary*)algo;      break;
    case kMichelID:       _alg_michel_id      = (BaseAlgMichelID*)algo;      break;
    case kMichelCluster:  _alg_michel_cluster = (BaseAlgMichelCluster*)algo; break;
    default:
      std::cerr << "\033[93m[ERROR]\033[00m "
		<< "Unidentified algorithm type: " << type << std::endl;
      throw MichelException();
    }

    _alg_v[type] = algo;
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
    // If nothing to be done, return
    if(_input_v.empty()) return;

    // make sure hit vectors provided & cluster list (another hit arrays) are consistent
    _used_hit_marker_v.resize(_all_hit_v.size(),false);
    for(size_t i=0; i<_used_hit_marker_v.size(); ++i) {
      if(_all_hit_v[i]._id != 2) 
	_used_hit_marker_v[i] = true;
      else
	_used_hit_marker_v[i] = false;
    }
    
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

	  std::cerr << "\033[93m[ERROR]\033[00m "
		    << "Found a hit index out of range! "
		    << hit_pt._id << " out of " << _used_hit_marker_v.size()
		    << std::endl;
	  throw MichelException();
	}
	_used_hit_marker_v[hit_pt._id] = true;
      }
    }

    //
    // Step 2 ... find muon/michel boundary
    //
    if(!_alg_boundary) throw MichelException();

    _watch.Start();
    for(auto& cluster : _output_v )

      // Find start point for michel
      cluster._boundary = _alg_boundary->Boundary(cluster);

    _alg_time_v [kBoundaryFinder] += _watch.RealTime();
    _alg_ctr_v  [kBoundaryFinder] += _output_v.size();

    //
    // Step 3 ... Identify michel/muon
    //
    if(!_alg_michel_id) throw MichelException();

    _watch.Start();
    for(auto& cluster : _output_v)

      cluster._michel = _alg_michel_id->Identify(cluster);

    _alg_time_v [kMichelID] += _watch.RealTime();
    _alg_ctr_v  [kMichelID] += _output_v.size();

    //
    // Step 4 ... Michel re-clustering
    //
    _watch.Start();
    if(_alg_michel_cluster) {      

      size_t ctr = 0;
      for(auto& cluster : _output_v) {

	std::vector<HitPt> available_hits_v;
	for( auto const& v : _used_hit_marker_v ) if(v) ++ctr;
	available_hits_v.reserve(ctr);
	for(size_t hit_index=0; hit_index<_all_hit_v.size(); ++hit_index)

	  if(!_used_hit_marker_v[hit_index]) available_hits_v.push_back(_all_hit_v[hit_index]);
	
	_alg_michel_cluster->Cluster(cluster._michel,available_hits_v);

      }
      
    }
    
    _alg_time_v [kMichelCluster] += _watch.RealTime();
    _alg_ctr_v  [kMichelCluster] += _output_v.size();
    
    //
    // Finally call analyze
    //
    for(auto& ana : _ana_v) ana->Analyze(_input_v,_output_v);
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

    double merge_time     = ( _alg_ctr_v[ kClusterMerger  ] ? _alg_time_v[ kClusterMerger  ] / ((double)(_alg_ctr_v[ kClusterMerger  ])) : 0 );
    double boundary_time  = ( _alg_ctr_v[ kBoundaryFinder ] ? _alg_time_v[ kBoundaryFinder ] / ((double)(_alg_ctr_v[ kBoundaryFinder ])) : 0 );
    double id_time        = ( _alg_ctr_v[ kMichelID       ] ? _alg_time_v[ kMichelID       ] / ((double)(_alg_ctr_v[ kMichelID       ])) : 0 );
    double recluster_time = ( _alg_ctr_v[ kMichelCluster  ] ? _alg_time_v[ kMichelCluster  ] / ((double)(_alg_ctr_v[ kMichelCluster  ])) : 0 );

    std::cout << std::endl
	      << "=================== Time Report =====================" << std::endl
	      << "  Merging        Algo Time: " << merge_time*1.e6     << " [us/cluster]" << std::endl
	      << "  Boudnary       Algo Time: " << boundary_time*1.e6  << " [us/cluster]" << std::endl
	      << "  Michel ID      Algo Time: " << id_time*1.e6        << " [us/cluster]" << std::endl
	      << "  Michel Cluster Algo Time: " << recluster_time*1.e6 << " [us/cluster]" << std::endl
	      << "=====================================================" << std::endl
	      << std::endl;
  }


}
#endif
