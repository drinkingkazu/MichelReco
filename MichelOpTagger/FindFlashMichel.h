/**
 * \file FindFlashMichel.h
 *
 * \ingroup MichelOpTagger
 * 
 * \brief Class def header for a class FindFlashMichel
 *
 * @author david caratelli
 */

/** \addtogroup MichelOpTagger

    @{*/

#ifndef LARLITE_FINDFLASHMICHEL_H
#define LARLITE_FINDFLASHMICHEL_H

#include <iostream>
#include "DataFormat/opflash.h"
#include "TTree.h"

namespace larlite {
  /**
     \class FindFlashMichel
   */
  class FindFlashMichel {
  
  public:

    /// Default constructor
    FindFlashMichel();

    /// Default destructor
    virtual ~FindFlashMichel(){}

    /**
     * @brief find candidate michel flash that matches a candidate decaying muon flash
     * @input vector<opflash> flashes -> all flashes in the event
     * @input opflash muon            -> flash associated w/ muon
     * @return position of flash for candidate michel in input vector. If none found returns -1
     */
    int FindMichelMatch(const std::vector<larlite::opflash>& flashes,
			const opflash& muon);

    /**
     * @brief get the compatibility for two flashes being muon-michel
     * @detail this function uses the distance between the two
     * as the metric to score a match (closer -> higher match)
     * @inpiut opflash muon   -> candidate muon flash
     * @inpiut opflash michel -> candidate michel flash
     * @return double -> score of flash compatibility
     */
    double MatchScoreDistance(const opflash& muon, const opflash& michel) const;

    /**
     * @brief get the compatibility for two flashes being muon-michel
     * @detail michel candidate w/ highest PE count has highest score
     * @inpiut opflash muon   -> candidate muon flash
     * @inpiut opflash michel -> candidate michel flash
     * @return double -> score of flash compatibility
     */
    double MatchScoreNPE(const opflash& michel) const;

    /**
     * @brief set time-window in which to search for michel
     * @input t -> time (usec)
     */
    void SetTimeWindow(double t) { _time_window = t; }

    /**
     * @brief set minimum num of PE for michel to be a match
     * @input pe -> photo electrons (PE)
     */
    void SetPEMin(double pe) { _PEmin = pe; }

    /**
     * @brief verobisty flag setter
     */
    void SetVerbose(bool on) { _verbose = on; }

    /**
     * @brief return tree filled by algorithm
     */
    TTree* GetTree() { return _tree; }

    /**
     * @brief initialize algorithm
     */
    void initialize();

  protected:

    /// verbosity flag
    bool _verbose;

    /// algorithm name
    std::string _name;

    /// search-time : how far in time search for the michel
    double _time_window;

    /// minimum number of PE allowed for match
    double _PEmin;

    TTree* _tree;
    int _n_compat;    // number of flashes within time allowed time-window
    double _max_PE;   // max. num of PE for all flahsed in the time-window
    double _min_d;    // min. distance for all flashes in the time-window;
    double _match_PE; // PE for matched flash (michel candidate)
    double _match_d;  // distance for matched flash (michel candidate)
    
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
