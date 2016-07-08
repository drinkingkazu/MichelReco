/**
 * \file NeutrinoFilter.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class NeutrinoFilter
 *
 * @author davidc1
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_NEUTRINOFILTER_H
#define LARLITE_NEUTRINOFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class NeutrinoFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class NeutrinoFilter : public ana_base{
  
  public:

    /// Default constructor
    NeutrinoFilter(){ _name="NeutrinoFilter"; _fout=0;}

    /// Default destructor
    virtual ~NeutrinoFilter(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

    int _tot;
    int _pass;
    
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
