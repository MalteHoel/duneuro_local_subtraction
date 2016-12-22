#ifndef DUNEURO_FITTED_TDCS_DRIVER_HH
#define DUNEURO_FITTED_TDCS_DRIVER_HH

#include <dune/common/parametertree.hh>

#include <duneuro/common/default_grids.hh>
#include <duneuro/common/dg_solver.hh>
#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/io/volume_conductor_reader.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/tes/fitted_tdcs_solver.hh>
#include <duneuro/tes/tdcs_driver_interface.hh>
#if HAVE_DUNE_SUBGRID
#include <duneuro/common/geometry_adaption.hh>
#endif

namespace duneuro
{
  template <FittedSolverType solverType, class VC, ElementType et, int degree>
  struct TDCSSelectFittedSolver;

  template <class VC, ElementType et, int degree>
  struct TDCSSelectFittedSolver<FittedSolverType::dg, VC, et, degree> {
    using SolverType = DGSolver<VC, et, degree>;
  };

  template <int d, ElementType elementType, bool geometryAdaption>
  class TDCSVolumeConductorStorage;

  template <int d, ElementType elementType>
  class TDCSVolumeConductorStorage<d, elementType, false>
  {
  public:
    using Type = VolumeConductor<typename DefaultGrid<d, elementType>::GridType>;

    explicit TDCSVolumeConductorStorage(const FittedDriverData<d>& data,
                                        const Dune::ParameterTree& config,
                                        DataTree dataTree = DataTree())
        : volumeConductor_(
              VolumeConductorReader<typename Type::GridType>::read(data, config, dataTree))
    {
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    std::shared_ptr<Type> volumeConductor_;
  };

#if HAVE_DUNE_SUBGRID
  // note: geometry adaption currently only available in 3d
  template <>
  class TDCSVolumeConductorStorage<3, ElementType::hexahedron, true>
  {
  public:
    using Type = VolumeConductor<typename GeometryAdaptedGrid<3>::GridType>;

    explicit TDCSVolumeConductorStorage(const FittedDriverData<3>& data,
                                        const Dune::ParameterTree& config,
                                        DataTree dataTree = DataTree())
        : adaptedGrid_(GeometryAdaptedGridReader<3>::read(config.sub("grid")))
        , volumeConductor_(make_geometry_adapted_volume_conductor<3>(
              std::move(adaptedGrid_.grid), std::move(adaptedGrid_.labels), config))
    {
    }

    std::shared_ptr<Type> get() const
    {
      assert(volumeConductor_);
      return volumeConductor_;
    }

  private:
    GeometryAdaptedGrid<3> adaptedGrid_;
    std::shared_ptr<Type> volumeConductor_;
  };
#endif

  template <int dim, ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption>
  struct FittedTDCSDriverTraits {
    using VCStorage = TDCSVolumeConductorStorage<dim, elementType, geometryAdaption>;
    using VC = typename VCStorage::Type;
    using Solver = typename TDCSSelectFittedSolver<solverType, VC, elementType, degree>::SolverType;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  };

  template <int dim, ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption = false>
  class FittedTDCSDriver : public TDCSDriverInterface<dim>
  {
  public:
    using Traits = FittedTDCSDriverTraits<dim, elementType, solverType, degree, geometryAdaption>;

    explicit FittedTDCSDriver(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : FittedTDCSDriver(FittedDriverData<dim>{}, config, dataTree)
    {
    }

    explicit FittedTDCSDriver(const FittedDriverData<dim>& data, const Dune::ParameterTree& config,
                              DataTree dataTree = DataTree())
        : config_(config)
        , volumeConductorStorage_(data, config.sub("volume_conductor"),
                                  dataTree.sub("volume_conductor"))
        , solver_(std::make_shared<typename Traits::Solver>(
              volumeConductorStorage_.get(),
              config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree()))
        , tdcsSolver_(volumeConductorStorage_.get(), solver_)
    {
    }

    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }

    virtual void solveTDCSForward(const PatchSet<double, dim>& patchSet, Function& solution,
                                  const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      tdcsSolver_.solve(patchSet, solution.cast<typename Traits::DomainDOFVector>(), config,
                        dataTree);
    }

    virtual void write(const Function& function, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        VTKWriter<typename Traits::VC, degree> writer(
            volumeConductorStorage_.get(), config.get<unsigned int>("subsampling", degree - 1));
        auto gradient_type = config.get<std::string>("gradient.type", "vertex");
        auto potential_type = config.get<std::string>("potential.type", "vertex");

        if (gradient_type == "vertex") {
          writer.addVertexDataGradient(
              tdcsSolver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "gradient_potential");

        } else {
          writer.addCellDataGradient(
              tdcsSolver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "gradient_potential");
        }
        if (potential_type == "vertex") {
          writer.addVertexData(tdcsSolver_, Dune::stackobject_to_shared_ptr(
                                                function.cast<typename Traits::DomainDOFVector>()),
                               "potential");

        } else {
          writer.addCellData(tdcsSolver_, Dune::stackobject_to_shared_ptr(
                                              function.cast<typename Traits::DomainDOFVector>()),
                             "potential");
        }
        writer.addCellData(std::make_shared<duneuro::TensorFunctor<typename Traits::VC>>(
            volumeConductorStorage_.get()));
        writer.write(config.get<std::string>("filename"), dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

  private:
    Dune::ParameterTree config_;
    typename Traits::VCStorage volumeConductorStorage_;
    std::shared_ptr<typename Traits::Solver> solver_;
    FittedTDCSSolver<typename Traits::Solver> tdcsSolver_;
  };
}

#endif // DUNEURO_FITTED_TDCS_DRIVER_HH
