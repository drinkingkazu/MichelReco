/**
 * \file Michel.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class Michel
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/
#ifndef MICHELCLUSTER_MICHEL_H
#define MICHELCLUSTER_MICHEL_H

#include "MichelTypes.h"
#include "HitPt.h"
#include <vector>
#include "ColorPrint.h"
namespace michel {

  /**
     \class Michel
  */
  class Michel : public ColorPrint,
		 public std::vector<michel::HitPt> {

  public:

    /// Default constructor
    Michel();
    
    Michel(const double charge,
	   const double energy,
	   const double length, 
	   const HitPt& start);
    
    /// Default destructor
    ~Michel(){}
    
    double _charge;
    double _energy;
    double _length;
    HitPt  _start;
    HitPt  _dir;

    std::string Diff(const Michel& rhs) const;
    
    void Dump() const;

    // store electron hit indices
    std::vector<size_t> _electron_hit_idx_v;
    // store photon hit clusters, each a list of hit indices
    std::vector< std::vector<size_t> > _photon_clus_v;
    
    
  };
}

#endif
/** @} */ // end of doxygen group 

