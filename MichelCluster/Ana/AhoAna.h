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

#include "TTree.h"

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
    
    TTree* _out_tree;
    double _largest_cluster_charge;
    double _n_hits_in_largest_cluster;
    double _n_hits_in_largest_cluster_michel;
    
    int _number_of_clusters;
    
    //Hits in the clusters
    std::vector<double> _Z;
    std::vector<double> _X;

    //Michels
    std::vector<double> michel_Z;
    std::vector<double> michel_X;
    
    //Etc...
    std::vector<double> charge_in_largest_cluster;
    std::vector<double> truncated_charge_in_largest_cluster;
    std::vector<double> truncated_dqds_in_largest_cluster;
    std::vector<double> s;
    
    
  };
}

#endif
/** @} */ // end of doxygen group 

