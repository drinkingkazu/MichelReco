/**
 * \file FilterMichels.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class FilterMichels
 *
 * @author davidc1
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_FILTERMICHELS_H
#define LARLITE_FILTERMICHELS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FilterMichels
     User custom analysis class made by SHELL_USER_NAME
   */
  class FilterMichels : public ana_base{
  
  public:

    /// Default constructor
    FilterMichels(){ _name="FilterMichels"; _fout=0; _producer="michel"; }

    /// Default destructor
    virtual ~FilterMichels(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetProducer(std::string s) { _producer = s; }

  protected:

    // require clusters by set producer name
    std::string _producer;

    int _totl_events;
    int _pass_events;
    
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
