/**
 * \file MichelRecoDriver.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class MichelRecoDriver
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/

#ifndef LARLITE_MICHELRECODRIVER_H
#define LARLITE_MICHELRECODRIVER_H

#include "Analysis/ana_base.h"
#include "MichelReco.h"
#include "MichelTypes.h"
#include <TTree.h>


//Backtracker
#include "MCComp/MCMatchAlg.h"



namespace larlite {
  /**
     \class MichelRecoDriver
   */
  class MichelRecoDriver : public ana_base{
  
  public:

    /// Default constructor
    MichelRecoDriver();

    /// Default destructor
    virtual ~MichelRecoDriver(){}

    /** IMPLEMENT in MichelRecoDriver.cc!
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MichelRecoDriver.cc! 
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MichelRecoDriver.cc! 
    */
    virtual bool finalize();

    /// Reco manager getter
    ::michel::MichelReco& Algo() { return _mgr; }

    /**
     * @brief Set boolean to decide if to save michels in output cluster
     */
    void saveMichelClusters(bool on) { _save_clusters = on; }
    
    /// Producer setter
    void SetClusterProducer(const std::string& name) { _producer = name; }

    /// Set electric field strength (kV/cm)
    void SetEField(double E) { _Efield = E; }

    /// set if to use MC info or not
    void SetUseMC(bool on) { _use_mc = on; }

  protected:

    /// Reco manager
    ::michel::MichelReco _mgr;

    /// Input cluster producer name string
    std::string _producer;

    // boolean. Use MC info?
    bool _use_mc;

    /// Output analysis TTree ptr
    TTree* _hit_tree;
    /// info for tree
    std::vector<double> _q_v;
    std::vector<double> _w_v;
    std::vector<double> _t_v;
    int _run;
    int _subrun;
    int _event;

    // mc tree
    TTree* _mc_tree;
    double _mc_energy;
    double _reco_energy;
    std::vector<double> _michel_hit_frac;
    
    /// boolean to select if to save Michel Clusters or not
    bool _save_clusters;

    /// electric field strength [ kV/cm]
    double _Efield;

  private:
    ::btutil::MCMatchAlg _BTAlg;


  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
