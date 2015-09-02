/**
 * \file CutOnMuonLength.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class CutOnMuonLength
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef CUTONMUONLENGTH_H
#define CUTONMUONLENGTH_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
  /**
     \class CutOnMuonLength
     User defined class CutOnMuonLength ... these comments are used to generate
     doxygen documentation!
  */
  class CutOnMuonLength : public BaseMichelAlgo {
    
  public:
    
    /// Default constructor
    CutOnMuonLength();
    
    /// Default destructor
    ~CutOnMuonLength(){};

    /// Event re-setter
    void EventReset(){};

    /**
     * @brief Use MichelCluster and surrounding hits to decide if this is really a michel
     * @input MichelCluster michel : the currently reconstructed michel object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    bool ProcessCluster(MichelCluster& michel,
			const std::vector<HitPt>& hits);

    /**
     *@brief set the minimum length for the muon
     */
    void SetMinMuonLength(double l) { _min_muon_length = l; }

  private:

    // minimum muon length. If the muon is shorter than this ignore this michel
    double _min_muon_length;
    
  };
}

#endif
/** @} */ // end of doxygen group 

