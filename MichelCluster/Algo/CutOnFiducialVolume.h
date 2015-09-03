/**
 * \file CutOnFiducialVolume.h
 *
 * \ingroup MichelCluster
 *
 * \brief Class def header for a class CutOnFiducialVolume
 *
 * @author david kaleko
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

    /// Setter to exclude a specific input wire range (in cm, cm units)
    void AddExcludedWireRange(std::pair<double, double> myrange)
    { _excluded_wire_ranges.push_back(myrange); }

    /// Setter to exclude a specific input time range (in cm, cm units)
    void AddExcludedTimeRange(std::pair<double, double> myrange)
    { _excluded_time_ranges.push_back(myrange); }

    /// Setter to exclude many wire ranges at once
    void SetExcludedWireRanges(std::vector<double> myranges_min,
			       std::vector<double> myranges_max);
      
      /// Setter to exclude many time ranges at once
      void SetExcludedTimeRanges(std::vector<double> myranges_min,
				 std::vector<double> myranges_max);
      
    /// Setter to show debug couts
    void SetDebug(bool flag) { _debug = flag; }

private:

    /// pairs of ([cm], [cm]) that indicate a (start, stop) list of dead wires
    /// code will exclude michels that have hits within _buffer [cm] of these boundaries
    /// (and of course no hits will exist in these exact ranges if wires are dead there)
    std::vector<std::pair<double, double>> _excluded_wire_ranges;
    std::vector<std::pair<double, double>> _excluded_time_ranges;

    /// buffer size. if buffer size is 3cm, and the range (12.3, 45.6) is excluded
    /// then if a michel has any hits in the range (12.3 - 3, 45.6 + 3) will be rejected
    double _wire_buffer_size;
    double _time_buffer_size;

    bool _debug;
};
}

#endif
/** @} */ // end of doxygen group

