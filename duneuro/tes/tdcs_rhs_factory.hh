#ifndef DUNEURO_TDCS_RHS_FACTORY_HH
#define DUNEURO_TDCS_RHS_FACTORY_HH

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <duneuro/tes/tdcs_point_rhs.hh>
#include <duneuro/tes/tdcs_rhs_interface.hh>

namespace duneuro
{
  struct FittedTdcsRHSFactory {
    template <class Vector, class Solver, class VC>
    static std::unique_ptr<
        TdcsRHSInterface<typename Solver::Traits::VolumeConductor::GridView, Vector>>
    create(const Solver& solver, const VC& vc, const Dune::ParameterTree& config)
    {
      auto type = config.get<std::string>("type", "point");
      if (type == "point") {
        return std::make_unique<FittedTdcsRHS<typename Solver::Traits::FunctionSpace::GFS, Vector>>(
            solver.functionSpace().getGFS());
      } else {
        DUNE_THROW(Dune::Exception, "unknown electrode type \"" << type << "\"");
      }
    }
  };

  struct UnfittedTdcsRHSFactory {
    template <class Vector, class Solver>
    static std::unique_ptr<
        TdcsRHSInterface<typename Solver::Traits::SubTriangulation::GridView, Vector>>
    create(const Solver& solver, typename Solver::Traits::SubTriangulation& subTriangulation,
           const Dune::ParameterTree& config)
    {
      auto type = config.get<std::string>("type", "point");
      if (type == "point") {
        return std::make_unique<UnfittedTdcsRHS<typename Solver::Traits::FunctionSpace::GFS, Vector,
                                                typename Solver::Traits::SubTriangulation>>(
            solver.functionSpace().getGFS(), subTriangulation,
            config.get<std::size_t>("compartments"));
      } else {
        DUNE_THROW(Dune::Exception, "unknown electrode type \"" << type << "\"");
      }
    }
  };
}

#endif // DUNEURO_TDCS_RHS_FACTORY_HH
