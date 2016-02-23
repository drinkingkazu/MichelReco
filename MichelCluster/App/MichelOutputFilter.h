/**
 * \file MichelOutputFilter.h
 *
 * \ingroup App
 * 
 * \brief Class def header for a class MichelOutputFilter
 *
 * @author david caratelli
 */

/** \addtogroup App

    @{*/

#ifndef LARLITE_MICHELOUTPUTFILTER_H
#define LARLITE_MICHELOUTPUTFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MichelOutputFilter
     \brief filter any event that does not have an output reconstructed michel
   */
  class MichelOutputFilter : public ana_base{
  
  public:

    /// Default constructor
    MichelOutputFilter();

    /// Default destructor
    virtual ~MichelOutputFilter(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /**
     * @brief set reconstructed michel cluster producer name
     */
    void setMichelClusterProducer(std::string s) { _michel_producer = s; }

  protected:

    /// producer name for michel clusters
    std::string _michel_producer;
    
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
