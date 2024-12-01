// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_LOGGED_TIMER_HH
#define DUNEURO_LOGGED_TIMER_HH

#include <dune/common/timer.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  class LoggedTimer
  {
  public:
    explicit LoggedTimer(DataTree dataTree) : dataTree_(dataTree)
    {
    }

    void lap(const std::string& name)
    {
      timer_.stop();
      dataTree_.set("time_" + name, timer_.lastElapsed());
      timer_.start();
    }

    void stop(const std::string& name = "")
    {
      timer_.stop();
      if (name == "") {
        dataTree_.set("time", timer_.elapsed());
      } else {
        dataTree_.set("time_" + name, timer_.elapsed());
      }
    }

  private:
    Dune::Timer timer_;
    DataTree dataTree_;
  };
}

#endif // DUNEURO_LOGGED_TIMER_HH
