#ifndef DUNEURO_DEPRECATED_HH
#define DUNEURO_DEPRECATED_HH

#include <string>

namespace duneuro
{
  inline void issueDeprecationWarning(const std::string& msg)
  {
    std::cout << "DEPRECATION WARNING:\n" << msg << std::endl;
  }
}

#endif // DUNEURO_DEPRECATED_HH
