/**
 * \file StoppingFilter.h
 *
 * \ingroup MichelElectron
 * 
 * \brief Class def header for a class StoppingFilter
 *
 * @author kathrynsutton
 */

/** \addtogroup MichelElectron

    @{*/

#ifndef LARLITE_STOPPINGFILTER_H
#define LARLITE_STOPPINGFILTER_H

#include "Analysis/ana_base.h"

//stolen from kas

namespace larlite {
  /**
     \class StoppingFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class StoppingFilter : public ana_base{
  
  public:

    /// Default constructor
    StoppingFilter();

    /// Default destructor
    virtual ~StoppingFilter(){}

    /** IMPLEMENT in StoppingFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in StoppingFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in StoppingFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    Int_t get_num_filtered() { return _kept_events; }

    void Stopping(bool on) { _stopping = on; }

    void Filter(bool on) { _filter = on; }

    void Pure(bool on) { _pure = on; }

    double GetEndEnergy() { return _E; }

    bool Passes() { return _pass; }

  protected:
    
    bool _stopping;

    bool _filter;
    
    bool _pass;
    
    int _kept_events;
    int _total_events;

    double _E;

    // boolean for a pure sample: i.e. ignore muons 
    // with end energy > 10 and < 100 MeV
    bool _pure;

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
