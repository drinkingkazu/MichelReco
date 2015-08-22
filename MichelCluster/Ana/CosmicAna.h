/**
 * \file SingleMuAna.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class SingleMuAna
 *
 * @author kazuhir0
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef COSMICANA_H
#define COSMICANA_H

#include "Fmwk/MichelAnaBase.h"

#include "Fmwk/MichelException.h"

#include "TTree.h"
#include "TH1F.h"
#include <math.h>

namespace michel {
  /**
     \class CosmicAna
     User defined class CosmicAna ... these comments are used to generate
     doxygen documentation!
  */
  class CosmicAna : public MichelAnaBase {
    
  public:
    
    /// Default constructor
    CosmicAna(){ _verbosity = msg::kNORMAL; }
    
    /// Default destructor
    ~CosmicAna(){}
    
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

    double _michel_clustered_charge;
    double _michel_n_hits;
    int    _number_of_clusters;

    int    _boundary;

    double _chi_at_boundary;
    double _mean_chi;
    double _rms_chi;
    double _lowest_chi;

    bool _has_michel;
    
    std::vector<double> _q_v;
    std::vector<double> _t_q_v;
    std::vector<double> _t_dqds_v;
    std::vector<double> _chi_v;
    
    std::vector<double> _s_v;
    
    std::vector<double> _Z;
    std::vector<double> _X;
    
    std::vector<double> _michel_Z;
    std::vector<double> _michel_X;

    // which part of the hit-list is the michel?
    int _forward;

    // also keep track of event information in the TTree
    int _run;
    int _subrun;
    int _event;
    int _clus_idx;
    
    
    double get_lowest(const std::vector<double>& data);
    double get_rms   (const std::vector<double>& data);
    double get_mean  (const std::vector<double>& data);
    void   clear_all();
    
    TH1F* _michel_hit_qs;
    
  };
}

#endif
/** @} */ // end of doxygen group 

