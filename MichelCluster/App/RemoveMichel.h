/**
 * \file RemoveMichel.h
 *
 * \ingroup MichelElectron
 * 
 * \brief Class def header for a class RemoveMichel
 *
 * @author kathrynsutton
 */

/** \addtogroup MichelElectron

    @{*/

#ifndef LARLITE_REMOVEMICHEL_H
#define LARLITE_REMOVEMICHEL_H

#include "Analysis/ana_base.h"

//stolen from kas

namespace larlite {
  /**
     \class RemoveMichel
     User custom analysis class made by SHELL_USER_NAME
  */
  class RemoveMichel : public ana_base{
  
  public:

    /// Default constructor
    RemoveMichel(){ _name="RemoveMichel"; _fout=0;}

    /// Default destructor
    virtual ~RemoveMichel(){}

    /** IMPLEMENT in RemoveMichel.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();
    
    /** IMPLEMENT in RemoveMichel.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RemoveMichel.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();
    
    Int_t get_num_filtered() { return kept_evts; }
    
  protected:
    
    size_t total_evts;
    size_t kept_evts;

    
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
