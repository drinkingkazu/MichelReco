/**
 * \file FindOpMichel.h
 *
 * \ingroup MichelOpTagger
 * 
 * \brief Class def header for a class FindOpMichel
 *
 * @author david caratelli
 */

/** \addtogroup MichelOpTagger

    @{*/

#ifndef LARLITE_FINDOPMICHEL_H
#define LARLITE_FINDOPMICHEL_H

#include <iostream>
#include "DataFormat/opflash.h"

namespace larlite {
  /**
     \class FindOpMichel
   */
  class FindOpMichel {
  
  public:

    /// Default constructor
    FindOpMichel();

    /// Default destructor
    virtual ~FindOpMichel(){}

    /**
     * @brief find candidate michel flash that matches a candidate decaying muon flash
     * @input vector<opflash> flashes -> all flashes in the event
     * @input opflash muon            -> flash associated w/ muon
     * @return position of flash for candidate michel in input vector. If none found returns -1
     */
    int FindMichelMatch(const std::vector<larlite::opflash>& flashes,
			const opflash& muon) const;

    /**
     * @brief get the compatibility for two flashes being muon-michel
     * @inpiut opflash muon   -> candidate muon flash
     * @inpiut opflash michel -> candidate michel flash
     * @return double -> score of flash compatibility
     */
    double MatchScore(const opflash& muon, const opflash& michel) const;

    /**
     * @brief set time-window in which to search for michel
     * @input t -> time (usec)
     */
    void SetTimeWindow(double t) { _time_window = t; }

  protected:

    /// algorithm name
    std::string _name;

    /// search-time : how far in time search for the michel
    double _time_window;
    
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
