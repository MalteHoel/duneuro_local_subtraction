// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNEURO_CUTFEM_GRIDOPERATOR_HH
#define DUNEURO_CUTFEM_GRIDOPERATOR_HH

#include <dune/istl/io.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/udg/pdelab/assembler/assembler.hh>
#include <dune/udg/pdelab/assembler/localassembler.hh>
#include <dune/udg/pdelab/ghostpenalty.hh>

namespace duneuro
{
  template <class GO, class ST, class ENP>
  class CutFEMGridOperator
  {
  public:
    using Traits = typename GO::Traits;

    //! Constructor for non trivial constraints
    explicit CutFEMGridOperator(const GO& go, std::shared_ptr<const ST> subTriangulation,
                                const ENP& edgeNormProvider, const Dune::ParameterTree& config)
        : go_(go)
        , subTriangulation_(subTriangulation)
        , edgeNormProvider_(edgeNormProvider)
        , ghost_penalty_(config.get<double>("ghost_penalty"))
        , conductivities_(config.get<std::vector<double>>("conductivities"))
    {
    }

    void fill_pattern(typename GO::Pattern& p) const
    {
      go_.fill_pattern(p);
    }

    void residual(const typename GO::Domain& x, typename GO::Range& r) const
    {
      go_.residual(x, r);
    }

    //! Assembler jacobian
    void jacobian(const typename GO::Domain& x, typename GO::Jacobian& a) const
    {
      go_.jacobian(x, a);
      Dune::UDG::add_ghost_penalty_cutfem(*subTriangulation_, go_.trialGridFunctionSpace(),
                                          edgeNormProvider_, ghost_penalty_, conductivities_, 0,
                                          Dune::PDELab::Backend::native(a));
    }

    const typename GO::Traits::TrialGridFunctionSpace& trialGridFunctionSpace() const
    {
      return go_.trialGridFunctionSpace();
    }

    const typename GO::Traits::TestGridFunctionSpace& testGridFunctionSpace() const
    {
      return go_.testGridFunctionSpace();
    }

    const typename GO::Traits::MatrixBackend& matrixBackend() const
    {
      return go_.matrixBackend();
    }

  private:
    const GO& go_;
    std::shared_ptr<const ST> subTriangulation_;
    ENP edgeNormProvider_;
    double ghost_penalty_;
    std::vector<double> conductivities_;
  };
}
#endif
