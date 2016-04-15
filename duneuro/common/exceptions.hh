#ifndef DUNEURO_EXCEPTIONS_HH
#define DUNEURO_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

namespace duneuro
{
  class Exception : public Dune::Exception
  {
  };

  class SourceModelException : public Exception
  {
  };
}

#endif // DUNEURO_EXCEPTIONS_HH
