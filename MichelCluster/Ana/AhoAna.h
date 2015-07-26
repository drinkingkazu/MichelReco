/**
 * \file AhoAna.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class AhoAna
 *
 * @author kazuhir0
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef AHOANA_H
#define AHOANA_H

#include "Fmwk/MichelAnaBase.h"

namespace michel {
  /**
     \class AhoAna
     User defined class AhoAna ... these comments are used to generate
     doxygen documentation!
  */
  class AhoAna : public MichelAnaBase {
    
  public:
    
    /// Default constructor
    AhoAna(){ _verbosity = msg::kNORMAL; }
    
    /// Default destructor
    ~AhoAna(){}
    
    /// Initialize
    void Initialize();
    
    /// Analyze
    void Analyze(const MichelClusterArray& input_cluster_v,
		 const MichelClusterArray& output_cluster_v);
    
    /// Event Reset
    void EventReset();
    
    /// Finalize
    void Finalize(TFile* fout);

    
  protected:

  };
}

#endif
/** @} */ // end of doxygen group 

