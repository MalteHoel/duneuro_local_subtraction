// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
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

  class IllegalArgumentException : public Exception
  {
  };
}

#endif // DUNEURO_EXCEPTIONS_HH
