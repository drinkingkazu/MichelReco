/**
 * \file EmptyMIDFilter.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class EmptyMIDFilter
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef EMPTYMIDFILTER_H
#define EMPTYMIDFILTER_H

#include "Fmwk/BaseMichelAlgo.h"

namespace michel {
/**
   \class EmptyMIDFilter
   User defined class EmptyMIDFilter ... these comments are used to generate
   doxygen documentation!
*/
class EmptyMIDFilter : public BaseMichelAlgo {

public:

    /// Default constructor
    EmptyMIDFilter();

    /// Default destructor
    ~EmptyMIDFilter() {};

    /// Event re-setter
    void EventReset() {};

    /**
     * @brief Use MichelCluster and surrounding hits to decide if this is really a michel
     * @input MichelCluster cluster : the currently reconstructed michel cluster object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    bool ProcessCluster(MichelCluster& cluster,
			const std::vector<HitPt>& hits);

};
}

#endif
/** @} */ // end of doxygen group

