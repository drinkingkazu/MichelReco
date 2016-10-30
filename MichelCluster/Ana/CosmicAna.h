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
    CosmicAna() { _verbosity = msg::kNORMAL; }

    /// Default destructor
    ~CosmicAna() {}

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

    // TTree where summary info from reconstructed Michels is stored for analysis
    TTree* _out_tree;

    double _michel_clustered_charge; ///< Total charge (in ADCs) for Michel electron
    int    _michel_n_hits;           ///< Number of hits tagged for Michel electron
    int    _electron_n_hits;         ///< Number of hits tagged for ionization-only segment
    int    _photon_n_hits;           ///< Number of hits tagged for photons associated to Michel
    int    _n_tagged_photons;        ///< Number of tagged photons associated to Michel
    int    _muon_n_hits;             ///< Number of hits associated to Muon cluster
    int    _boundary;                ///< Hit index for boundary between muon and Michel

    std::vector<double> _q_v;
    // variables indicating linearity value (chi^2) calculated for various parameters
    double _chi_at_boundary;         ///< Measured linearity at Muon-Michel boundary
    double _mean_chi;                ///< Mean linearity for all hits
    double _rms_chi;                 ///< RMS of linearity
    double _lowest_chi;              ///< minimum linearity value
    double _mean_chi_michel;         ///< mean linearity for michel hits
    double _mean_chi_muon;           ///< mean linearity for muon hits

    bool _has_michel;

    // list of indices of hits in Michel cluster associated with:
    // the various photon clusters
    std::vector<int> _photon_clus_v;
    // the electron cluster
    std::vector<int> _electron_clus;
    
    // vector profiles of hits in muon-michel
    std::vector<double> _t_q_v;      ///< Truncated mean charge profile
    std::vector<double> _t_dqds_v;   ///< Truncated dQ/ds charge profile
    std::vector<double> _chi_v;      ///< hit linearity profile
    std::vector<double> _dirs_v;     ///< vector of 2D directions between consecutive hits
    std::vector<double> _s_v;        ///< displacement along trajectory vector

    std::vector<double> _Z;          ///< vector of all Z positions
    std::vector<double> _X;          ///< vector of all X positions
    std::vector<int>    _idx;        ///< vector of hit indices (references Gaus hit his)

    std::vector<double> _michel_Z;   ///<
    std::vector<double> _michel_X;
    std::vector<double> _michel_Q;

    // which part of the hit-list is the michel?
    int _forward;

    // also keep track of event information in the TTree
    int _run;
    int _subrun;
    int _event;
    int _clus_idx;
    
    double _michel_start_Z, _michel_start_X;

    double _michel_dir_x, _michel_dir_z;

    // mean Q along muon truncating the averaging
    double _mean_q_muon;
    

    std::vector<double> _slope_v;


    double get_lowest(const std::vector<double>& data);
    double get_rms   (const std::vector<double>& data,
		      const double& avg);
    double get_mean  (const std::vector<double>& data);
    double get_mean_between_bounds(const std::vector<double>& data,
				   const double& minVal,
				   const double& maxVal);
    void   clear_all();

    TH1F* _michel_hit_qs;

};
}

#endif
/** @} */ // end of doxygen group

