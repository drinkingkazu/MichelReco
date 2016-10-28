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
#include "Fmwk/MichelRecoManager.h"
#include "Fmwk/MichelTypes.h"
#include <TStopwatch.h>
#include <TTree.h>
#include <fstream>


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
    ::michel::MichelRecoManager& GetManager() { return _mgr; }

    /**
     * @brief Set boolean to decide if to save michels in output cluster
     */
    void saveMichelClusters(bool on) { _save_clusters = on; }
    
    /// Producer setter
    void SetClusterProducer(const std::string& name) { _producer = name; }

    /// Specifier for a specific plane reconstruction
    void SetPlane(unsigned int p) { _reco_plane.resize(p+1,false); _reco_plane[p]=true; }

    /// Set the minimum number of hits for a cluster to
    /// be considered as a candidate MichelCluster
    void SetMinClusSize(int n) { _minClusSize = n; }

    /// filter events
    void FilterEvents(bool on) { _filter_events = on; }

    /// set channel-by-channel gain corrections
    void SetNChannels(int n) { _chchgain.resize( n ); }
    /// set a specific channel's gain
    void SetChGain(const int& ch, const double& g);
    /// set average gain so that charge is preserved in ADCs x Ticks
    void SetAvgGain(double g) { _avg_gain = g; }

    // set output txt file
    void set_outtxt_file(std::string s) { _out_txtfile = s; }

  protected:

    /// Stopwatch for event time profiling
    TStopwatch _event_watch; ///< For profiling
    double _event_time;
    size_t _event_ctr;
    
    /// Reco manager
    ::michel::MichelRecoManager _mgr;

    /// Input cluster producer name string
    std::string _producer;

    /// Option to set specific-plane-only reco
    std::vector<bool> _reco_plane;

    /// Output analysis TTree ptr
    TTree* _hit_tree;
    TTree* _evt_tree;
    /// info for tree
    std::vector<double> _q_v;
    std::vector<double> _w_v;
    std::vector<double> _t_v;
    std::vector<double> _p_v;
    
    int _run;
    int _subrun;
    int _event;

    std::vector<michel::EventID> _events_info;

    /// filter events
    bool _filter_events;

    /// boolean to select if to save Michel Clusters or not
    bool _save_clusters;

    /// minimum cluster size for an input cluster to be considered
    int _minClusSize;

    // channel-by-channel gain corrections
    std::vector<double> _chchgain;
    // average gain
    double _avg_gain;

    ofstream outfile;
    std::string _out_txtfile;

  private:

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
