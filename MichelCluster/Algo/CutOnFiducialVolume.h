/**
 * \file CutOnFiducialVolume.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class CutOnFiducialVolume
 *
 * @author david caratelli
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHEL_CUTONFIDUCIALVOLUME_H
#define MICHEL_CUTONFIDUCIALVOLUME_H

#include "Fmwk/BaseMichelAlgo.h"
#include "math.h"

namespace michel {
/**
   \class CutOnFiducialVolume
   User defined class CutOnFiducialVolume ... these comments are used to generate
   doxygen documentation!
*/
class CutOnFiducialVolume : public BaseMichelAlgo {

public:

    /// Default constructor
    CutOnFiducialVolume();

    /// Default destructor
    ~CutOnFiducialVolume() {};

    /// Event re-setter
    void EventReset() {};

    /**
     * @brief Use MichelCluster and surrounding hits to decide if this is really a michel
     * @input MichelCluster michel : the currently reconstructed michel object
     * @input std::vector<HitPt> hits : all hits in the event
     * @return boolean : is this truly a michel or not
     */
    bool ProcessCluster(MichelCluster& michel,
                        const std::vector<HitPt>& hits);

private:

};
}

#endif
/** @} */ // end of doxygen group

