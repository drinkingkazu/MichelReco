/**
 * \file MichelException.h
 *
 * \ingroup Michel
 * 
 * \brief Class def header for a class MichelException
 *
 * @author kazuhiro
 */

/** \addtogroup Michel

    @{*/
#ifndef MICHELCLUSTEREXCEPTION_H
#define MICHELCLUSTEREXCEPTION_H

#include <iostream>

namespace michel {
  /**
     \class MichelException
  */
  class MichelException : public std::exception {
    
  public:
    
    /// Default constructor
    MichelException(){}
    
    /// Default destructor
    ~MichelException(){}

    
    
  };
}

#endif
/** @} */ // end of doxygen group 

