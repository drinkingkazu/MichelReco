/**
 * \file MichelClusterUnitTest.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class MichelClusterUnitTest
 *
 * @author davidkaleko
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_MICHELCLUSTERUNITTEST_H
#define LARLITE_MICHELCLUSTERUNITTEST_H

#include "Analysis/ana_base.h"
#include "DataFormat/cluster.h"
 
namespace larlite {
  /**
     \class MichelClusterUnitTest
     User custom analysis class made by SHELL_USER_NAME
   */
  class MichelClusterUnitTest : public ana_base{
  
  public:

    /// Default constructor
    MichelClusterUnitTest(){ _name="MichelClusterUnitTest"; _fout=0;}

    /// Default destructor
    virtual ~MichelClusterUnitTest(){}

    /** IMPLEMENT in MichelClusterUnitTest.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MichelClusterUnitTest.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MichelClusterUnitTest.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    
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
