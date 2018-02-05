#include <duneuro/common/deprecated.hh>

#include <iostream>

namespace duneuro
{
  void issueDeprecationWarning(const std::string& msg)
  {
    std::cout << "DEPRECATION WARNING:\n" << msg << std::endl;
  }
}
