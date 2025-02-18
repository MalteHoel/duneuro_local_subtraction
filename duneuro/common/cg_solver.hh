// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_CG_SOLVER_HH
#define DUNEURO_CG_SOLVER_HH

#include <dune/common/math.hh>
#include <dune/common/version.hh>
#ifndef DUNE_VERSION_NEWER
#define DUNE_VERSION_NEWER(a,b,c) DUNE_VERSION_GTE(a,b,c)
#endif

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

#include <duneuro/common/assembler.hh>
#include <duneuro/common/convection_diffusion_cg_default_parameter.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/kdtree.hh>
#include <duneuro/common/linear_problem_solver.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/random.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class VolumeConductor, class BCType, unsigned int degree, ElementType elementType>
  struct CGFunctionSpaceTraits;

  template <class VolumeConductor, class BCType, unsigned int degree>
  struct CGFunctionSpaceTraits<VolumeConductor, BCType, degree, ElementType::tetrahedron> {
    using Type =
        Dune::PDELab::CGSpace<typename VolumeConductor::GridType, typename VolumeConductor::ctype,
                              degree, BCType, Dune::GeometryType::simplex,
                              Dune::PDELab::MeshType::conforming, Dune::SolverCategory::sequential>;
  };

  template <class VolumeConductor, class BCType, unsigned int degree>
  struct CGFunctionSpaceTraits<VolumeConductor, BCType, degree, ElementType::hexahedron> {
    using Type =
        Dune::PDELab::CGSpace<typename VolumeConductor::GridType, typename VolumeConductor::ctype,
                              degree, BCType, Dune::GeometryType::cube,
                              Dune::PDELab::MeshType::conforming, Dune::SolverCategory::sequential>;
  };

  template <class VC, ElementType elementType, unsigned int degree, class DF = double,
            class RF = double, class JF = double>
  struct CGSolverTraits {
    static const int dimension = VC::dim;
    static const bool isFitted = true;
    using VolumeConductor = VC;
    using GridView = typename VC::GridView;
    using CoordinateFieldType = typename VC::ctype;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using Problem = ConvectionDiffusionCGDefaultParameter<VC>;
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    using BoundaryCondition = Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem>;
    using FunctionSpace =
        typename CGFunctionSpaceTraits<VC, BoundaryCondition, degree, elementType>::Type;
    using DomainDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, DF>;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<typename FunctionSpace::GFS, RF>;
    using LocalOperator =
        Dune::PDELab::ConvectionDiffusionFEM<Problem, typename FunctionSpace::FEM>;
    using Assembler = GalerkinGlobalAssembler<FunctionSpace, LocalOperator, DF, RF, JF>;
    using LinearSolver =
        LinearProblemSolver<typename Assembler::GO, DomainDOFVector, RangeDOFVector>;
  };

  template <class VC, ElementType elementType, unsigned int degree, class DF = double,
            class RF = double, class JF = double>
  class CGSolver
  {
  public:
    using Traits = CGSolverTraits<VC, elementType, degree, DF, RF, JF>;

    CGSolver(std::shared_ptr<const VC> volumeConductor,
             std::shared_ptr<const typename Traits::ElementSearch> search,
             const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : volumeConductor_(volumeConductor)
        , search_(search)
        , problem_(volumeConductor)
        , dirichletExtension_(volumeConductor->gridView(), problem_)
        , boundaryCondition_(volumeConductor->gridView(), problem_)
        , functionSpace_(const_cast<typename VC::GridType&>(volumeConductor->grid()),
                         boundaryCondition_)
        // const_cast due to misused non-const reference in pdelab boilerplate
        , localOperator_(problem_, config.get<unsigned int>("intorderadd", 0))
        , assembler_(functionSpace_, localOperator_, elementType == ElementType::hexahedron ?
                                                         (1 << VC::dim) + 1 :
#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS, 2, 10)
                                                         Dune::power(3, int(VC::dim)))
#else
                                                         Dune::StaticPower<3, VC::dim>::power)
#endif
        , linearSolver_(assembler_.getGO(), config)
    {
      dataTree.set("degree", degree);
      dataTree.set("element_type", to_string(elementType));
    }

    template <typename SolverBackend>
    void solve(SolverBackend& solverBackend, const typename Traits::RangeDOFVector& rightHandSide,
               typename Traits::DomainDOFVector& solution, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      randomize_uniform(Dune::PDELab::Backend::native(solution), DF(-1.0), DF(1.0));
      linearSolver_.apply(solverBackend, solution, rightHandSide, config, dataTree);
      dataTree.set("time", timer.elapsed());
    }

    const typename Traits::FunctionSpace& functionSpace() const
    {
      return functionSpace_;
    }

    std::shared_ptr<const typename Traits::VolumeConductor> volumeConductor() const
    {
      return volumeConductor_;
    }

    std::shared_ptr<const typename Traits::ElementSearch> elementSearch() const
    {
      return search_;
    }

  private:
    std::shared_ptr<const typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<const typename Traits::ElementSearch> search_;
    typename Traits::Problem problem_;
    typename Traits::DirichletExtension dirichletExtension_;
    typename Traits::BoundaryCondition boundaryCondition_;
    typename Traits::FunctionSpace functionSpace_;
    typename Traits::LocalOperator localOperator_;
    typename Traits::Assembler assembler_;
    typename Traits::LinearSolver linearSolver_;

    template <class V>
    friend struct MakeDOFVectorHelper;
  };
}

#endif // DUNEURO_CG_SOLVER_HH
