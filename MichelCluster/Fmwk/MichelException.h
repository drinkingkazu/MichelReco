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
#include <exception>
namespace michel {
  /**
     \class MichelException
  */
  class MichelException : public std::exception {
    
  public:

    MichelException(const std::string& msg="");
    
    virtual ~MichelException() throw(){}

    virtual const char* what() const throw();

  private:

    std::string _msg;
    
  };
}

#endif
/** @} */ // end of doxygen group 

