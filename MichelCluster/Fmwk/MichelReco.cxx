#ifndef MICHELRECO_CXX
#define MICHELRECO_CXX

#include "MichelReco.h"
#include "MichelException.h"
namespace michel {

  //-----------------------------------------------------------------
  MichelReco::MichelReco()
  //-----------------------------------------------------------------
    : _verbosity  ( msg::kNORMAL )
    , _alg_v      ( kAlgoTypeMax, nullptr )
    , _alg_time_v ( kAlgoTypeMax, 0.      )
    , _alg_ctr_v  ( kAlgoTypeMax, 0       )
  {}
  
  //-----------------------------------------------------------------
  void MichelReco::Append(const std::vector<michel::HitPt>& hit_v)
  //-----------------------------------------------------------------
  {
    if(hit_v.size() < _min_nhits) return;
    MichelCluster cluster(_min_nhits, _d_cutoff);
    cluster.SetVerbosity(_verbosity);
    cluster.SetHits(hit_v);
    _input_v.emplace_back(cluster);
  }
  
  //-----------------------------------------------------------------
  void MichelReco::SetAlgo(const AlgoType_t type,
			   BaseMichelAlgo* algo)
  //-----------------------------------------------------------------
  {
    switch(type) {
    case kClusterMerger:  _alg_merge    = (BaseAlgMerger*)algo;     break;
    case kBoundaryFinder: _alg_boundary = (BaseAlgBoundary*)algo;   break;
    case kMichelFinder:   _alg_michel   = (BaseAlgIdentifier*)algo; break;
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
    if(!_alg_michel) throw MichelException();

    _watch.Start();
    for(auto& cluster : _output_v)

      cluster._michel = _alg_michel->Identify(cluster);

    _alg_time_v [kMichelFinder] += _watch.RealTime();
    _alg_ctr_v  [kMichelFinder] += _output_v.size();
    
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
  }

  //-----------------------------------------------------------------
  void MichelReco::Finalize(TFile *fout)
  //-----------------------------------------------------------------
  {
    for(auto& ana : _ana_v) ana->Finalize(fout);

    std::cout << "=================== Time Report =====================" << std::endl
	      << "Merging    Algo Time: " << _alg_time_v[ kClusterMerger  ] / ((double)(_alg_ctr_v[kClusterMerger]))  << " [s/event]" << std::endl
	      << "Boudnary   Algo Time: " << _alg_time_v[ kBoundaryFinder ] / ((double)(_alg_ctr_v[kBoundaryFinder])) << " [s/event]" << std::endl
	      << "Identifier Algo Time: " << _alg_time_v[ kMichelFinder   ] / ((double)(_alg_ctr_v[kMichelFinder]))   << " [s/event]" << std::endl
	      << "=====================================================" << std::endl
	      << std::endl;
  }


}
#endif
