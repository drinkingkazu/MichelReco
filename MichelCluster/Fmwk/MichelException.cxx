#ifndef MICHELCLUSTEREXCEPTION_CXX
#define MICHELCLUSTEREXCEPTION_CXX

#include "MichelException.h"
#include "MichelTypes.h"
namespace michel {

  MichelException::MichelException(const std::string& msg)
    : std::exception()
  {
    _msg += msg::kColorPrefix[msg::kEXCEPTION].c_str();
    _msg += msg::kStringPrefix[msg::kEXCEPTION].c_str();
    _msg += "\033[0m";
    _msg += msg;
    _msg += "\n";
  }

  const char* MichelException::what() const throw()
  { return _msg.c_str(); }


}

#endif
