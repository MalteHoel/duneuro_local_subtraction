#ifndef DUNEURO_UDG_TDCS_SOLVER_HH
#define DUNEURO_UDG_TDCS_SOLVER_HH

#include <duneuro/common/udg_solver.hh>
#include <duneuro/tes/tdcs_patch_udg_assembler.hh>

namespace duneuro
{
  template <class ST, int compartments, int degree, class DF, class RF, class JF>
  struct UDGTDCSSolverTraits {
    static const unsigned int dimension = ST::dim;
    using Solver = UDGSolver<ST, compartments, degree, TDCSPatchUDGParameter<typename ST::GridView>,
                             DF, RF, JF>;
    using SubTriangulation = ST;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using CoordinateFieldType = typename ST::ctype;
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGTDCSSolver
  {
  public:
    using Traits = UDGTDCSSolverTraits<ST, compartments, degree, DF, RF, JF>;

    UDGTDCSSolver(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                  std::shared_ptr<typename Traits::Solver> solver,
                  const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , solver_(solver)
        , rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
        , config_(config)
    {
    }

    void solve(typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      // assemble right hand side
      Dune::Timer timer;
      *rightHandSideVector_ = 0.0;
      TDCSPatchUDGAssembler<typename Traits::SubTriangulation, typename Traits::FunctionSpace,
                            typename Traits::RangeDOFVector>
          assembler(solver_->problem(), subTriangulation_, solver_->functionSpace(), config_);
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

    const typename Traits::SubTriangulation& subTriangulation() const
    {
      return *subTriangulation_;
    }

  private:
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::unique_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;

    Dune::ParameterTree config_;
  };
}

#endif // DUNEURO_UDG_TDCS_SOLVER_HH
