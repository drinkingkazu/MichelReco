/**
 * \file MichelAnaBase.h
 *
 * \ingroup MichelCluster
 * 
 * \brief Class def header for a class MichelAnaBase
 *
 * @author kazuhiro
 */

/** \addtogroup MichelCluster

    @{*/
#ifndef MICHELANABASE_H
#define MICHELANABASE_H

#include "MichelCluster.h"
#include <TFile.h>

namespace michel {
  /**
     \class MichelAnaBase
     User defined class MichelAnaBase ... these comments are used to generate
     doxygen documentation!
  */
  class MichelAnaBase{
    
  public:
    
    /// Default constructor
    MichelAnaBase(){ _verbosity = msg::kNORMAL; }
    
    /// Default destructor
    virtual ~MichelAnaBase(){}

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level)
    { _verbosity = level; }
    
    /// Initialize
    virtual void Initialize() = 0;

    /// Analyze
    virtual void Analyze(const MichelClusterArray& input_cluster_v,
			 const MichelClusterArray& output_cluster_v) = 0;
    /// Event Reset
    virtual void EventReset() = 0;

    /// Finalize
    virtual void Finalize(TFile* fout) = 0;

  protected:

    /// Verbosity level
    msg::MSGLevel_t _verbosity;
    
  };
}

#endif
/** @} */ // end of doxygen group 

