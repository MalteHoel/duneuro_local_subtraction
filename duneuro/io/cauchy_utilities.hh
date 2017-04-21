#ifndef DUNEURO_CAUCHY_COMMON_HH
#define DUNEURO_CAUCHY_COMMON_HH

#include <iostream>
#include <limits>

#include <dune/common/exceptions.hh>

namespace duneuro
{
  namespace CauchyDetail
  {
    template <char c>
    struct AllEqual {
      bool operator()(char o)
      {
        return c == o;
      }
    };

    static inline void strip(std::string& str)
    {
      str.erase(str.find_last_not_of(" \t") + 1);
      str.erase(0, str.find_first_not_of(" \t"));
    }

    static inline void ignoreSeparatorLine(std::istream& stream, unsigned int count = 1)
    {
      for (unsigned int i = 0; i < count; ++i) {
        std::string str;
        std::getline(stream, str, '\n');
        std::size_t index = str.find_first_not_of('=');
        if (index != std::string::npos) {
          DUNE_THROW(Dune::IOError, "expected separator line but found char: " << str[index]);
        }
      }
    }

    static inline void readLine(std::istream& stream, const std::string& str)
    {
      std::string line;
      std::getline(stream, line, '\n');
      if (line != str) {
        DUNE_THROW(Dune::IOError, "expected \"" << str << "\" but found \"" << line << "\"");
      }
    }

    static inline void skipWhitespaces(std::istream& stream)
    {
      stream.ignore(std::numeric_limits<std::size_t>::max(), ' ');
      stream.ignore(std::numeric_limits<std::size_t>::max(), '\n');
      stream.ignore(std::numeric_limits<std::size_t>::max(), ' ');
    }
  }
}
#endif // DUNEURO_CAUCHY_COMMON_HH
