// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later OR GPL-3.0-or-later
#ifndef DUNEURO_ASSEMBLER_HH
#define DUNEURO_ASSEMBLER_HH

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

namespace duneuro
{
  /**
   * Similar to Dune::PDELab::GalerkinGlobalAssembler, but let the user choose
   * the field types of the domain, range and jacobian
   */
  template <typename FS, typename LOP, typename DF = double, typename RF = double,
            typename JF = double>
  class GalerkinGlobalAssembler
  {
  public:
    // export types
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    using GO = Dune::PDELab::GridOperator<typename FS::GFS, typename FS::GFS, LOP, MBE, DF, RF, JF,
                                          typename FS::CC, typename FS::CC>;
    using MAT = typename GO::Jacobian;

    GalerkinGlobalAssembler(const FS& fs, LOP& lop, const std::size_t nonzeros)
    {
      gop = std::make_shared<GO>(fs.getGFS(), fs.getCC(), fs.getGFS(), fs.getCC(), lop,
                                 MBE(nonzeros));
    }

    // return grid reference
    GO& getGO()
    {
      return *gop;
    }

    // return grid reference const version
    const GO& getGO() const
    {
      return *gop;
    }

    GO& operator*()
    {
      return *gop;
    }

    GO* operator->()
    {
      return gop.operator->();
    }

    const GO& operator*() const
    {
      return *gop;
    }

    const GO* operator->() const
    {
      return gop.operator->();
    }

  private:
    std::shared_ptr<GO> gop;
  };
}

#endif // DUNEURO_ASSEMBLER_HH
