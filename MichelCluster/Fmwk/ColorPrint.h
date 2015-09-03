/**
 * \file ColorPrint.h
 *
 * \ingroup Fmwk
 * 
 * \brief Message utility class definition header file.
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup Fmwk

    @{*/


#ifndef MICHELCLUSTER_COLORPRINT_H
#define MICHELCLUSTER_COLORPRINT_H

#include <cstdio>
#include <iostream>
#include "MichelTypes.h"

namespace michel {
  /**
     \class Message
     \brief Utility class used to show formatted message on the screen.
  */
  class ColorPrint {
    
  public:
    
    /// Default constructor
    ColorPrint() : _verbosity(msg::kNORMAL) {}
    
    /// Default destructor
    virtual ~ColorPrint(){}

    /// Verbosity setter
    void SetVerbosity(msg::MSGLevel_t level) { _verbosity = level; }

    /// Print function
    void Print(msg::MSGLevel_t level, const std::string& header, const std::string& msg) const;

    /// Print function
    void Print(msg::MSGLevel_t level, const std::string& msg) const;

  protected:
    
    /// Verbosity
    msg::MSGLevel_t _verbosity;
    
  };
}
#endif
  
/** @} */ // end of doxygen group ColorPrint
