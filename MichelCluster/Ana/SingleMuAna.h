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
#ifndef SINGLEMUANA_H
#define SINGLEMUANA_H

#include "Fmwk/MichelAnaBase.h"

#include "Fmwk/MichelException.h"

#include "TTree.h"

namespace michel {
  /**
     \class SingleMuAna
     User defined class SingleMuAna ... these comments are used to generate
     doxygen documentation!
  */
  class SingleMuAna : public MichelAnaBase {
    
  public:
    
    /// Default constructor
    SingleMuAna(){ _verbosity = msg::kNORMAL; }
    
    /// Default destructor
    ~SingleMuAna(){}
    
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
    std::vector<double> covariance_in_largest_cluster;

    double chi_at_boundary;
    double mean_chi;
    double lowest_chi;
    double rms_chi;


    
    
    int boundary;
    int IMAX;
    int IMIN;


    double get_lowest(const std::vector<double>& data);
    double get_rms(const std::vector<double>& data);
    double get_mean(const std::vector<double>& data);

    
  };
}

#endif
/** @} */ // end of doxygen group 

