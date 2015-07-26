#ifndef AHOANA_CXX
#define AHOANA_CXX

#include "AhoAna.h"

namespace michel {

  /// Initialize
  void AhoAna::Initialize()
  {
    std::cout << "Initialize\n";
    
  }
  
  /// Analyze
  void AhoAna::Analyze(const MichelClusterArray& input_cluster_v,
		       const MichelClusterArray& output_cluster_v)
  {
    for(auto& output : output_cluster_v) {
      output._michel.Dump();
    }
    
  }
  
  /// Event Reset
  void AhoAna::EventReset()
  {
    
  }
  
  /// Finalize
  void AhoAna::Finalize(TFile* fout)
  {
    
  }
  
}
#endif

