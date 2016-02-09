/**
 * \file StoppingMuons.h
 *
 * \ingroup MichelElectron
 * 
 * \brief Class def header for a class StoppingMuons
 *
 * @author kathrynsutton
 */

/** \addtogroup MichelElectron

    @{*/

#ifndef LARLITE_STOPPINGMUONS_H
#define LARLITE_STOPPINGMUONS_H

#include "Analysis/ana_base.h"

//stolen from kas

namespace larlite {
  /**
     \class StoppingMuons
     User custom analysis class made by SHELL_USER_NAME
   */
  class StoppingMuons : public ana_base{
  
  public:

    /// Default constructor
    StoppingMuons(){ _name="StoppingMuons"; _fout=0; _flip=false;}

    /// Default destructor
    virtual ~StoppingMuons(){}

    /** IMPLEMENT in StoppingMuons.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in StoppingMuons.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in StoppingMuons.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    Int_t get_num_filtered() { return kept_evts; }

    void flip_filter(bool doit) { _flip = doit; }

  protected:
    
    //TTree* n_events_tree;
    //Whether to flip the output of the filter
    bool _flip;

    size_t total_evts;
    size_t kept_evts;
    
    TTree* michel_filter_tree;
    
    double _michel_energy;
    double _michel_charge;
    double _michel_x;
    double _michel_det;
    double _lifetime_correction;

    TTree* event_filter_tree;
    
    int _run,_subrun,_event;
    int _n_michel;
    
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
