#ifndef DUNEURO_FITTED_TDCS_SOLVER_HH
#define DUNEURO_FITTED_TDCS_SOLVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/tes/patch_set.hh>
#include <duneuro/tes/tdcs_patch_dg_assembler.hh>

namespace duneuro
{
  template <class S>
  struct FittedTDCSSolverTraits {
    static const unsigned int dimension = S::Traits::dimension;
    using Solver = S;
    using VolumeConductor = typename S::Traits::VolumeConductor;
    using FunctionSpace = typename S::Traits::FunctionSpace;
    using DomainDOFVector = typename S::Traits::DomainDOFVector;
    using RangeDOFVector = typename S::Traits::RangeDOFVector;
    using CoordinateFieldType = typename VolumeConductor::ctype;
  };

  template <class S>
  class FittedTDCSSolver
  {
  public:
    using Traits = FittedTDCSSolverTraits<S>;

    FittedTDCSSolver(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                     std::shared_ptr<typename Traits::Solver> solver)
        : volumeConductor_(volumeConductor)
        , solver_(solver)
        , rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
    {
    }

    void solve(const PatchSet<typename Traits::CoordinateFieldType, Traits::VolumeConductor::dim>&
                   patchSet,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      TDCSPatchDGAssembler<typename Traits::VolumeConductor, typename Traits::FunctionSpace,
                           typename Traits::RangeDOFVector>
          assembler(patchSet, volumeConductor_, solver_->functionSpace(), config);
      assembler.assembleRightHandSide(*rightHandSideVector_);
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();
      // solve system
      solver_->solve(*rightHandSideVector_, solution, config.sub("solver"),
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solve", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return solver_->functionSpace();
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_FITTED_TDCS_SOLVER_HH
