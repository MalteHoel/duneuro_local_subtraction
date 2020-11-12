#ifndef DUNEURO_FITTED_TDCS_DRIVER_HH
#define DUNEURO_FITTED_TDCS_DRIVER_HH

#include <memory>

#include <dune/common/parametertree.hh>

#include <duneuro/common/default_grids.hh>
#include <duneuro/common/dg_solver.hh>
#include <duneuro/common/dg_solver_backend.hh>
#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/volume_conductor_storage.hh>
#include <duneuro/io/fitted_tensor_vtk_functor.hh>
#include <duneuro/io/volume_conductor_reader.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/tes/tdcs_driver_interface.hh>
#include <duneuro/tes/tdcs_patch_dg_parameter.hh>
#if HAVE_DUNE_SUBGRID
#include <duneuro/common/geometry_adaption.hh>
#endif

namespace duneuro
{
  template <FittedSolverType solverType, class VC, ElementType et, int degree>
  struct TDCSSelectFittedSolver;

  template <class VC, ElementType et, int degree>
  struct TDCSSelectFittedSolver<FittedSolverType::dg, VC, et, degree> {
    using ProblemType = TDCSPatchDGParameter<VC>;
    using SolverType = DGSolver<VC, et, degree, ProblemType>;
    using SolverBackendType = DGSolverBackend<SolverType, et>;
  };

  template <int dim, ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption>
  struct FittedTDCSDriverTraits {
    using VCStorage = VolumeConductorStorage<dim, elementType, geometryAdaption>;
    using VC = typename VCStorage::Type;
    using Solver = typename TDCSSelectFittedSolver<solverType, VC, elementType, degree>::SolverType;
    using SolverBackend =
        typename TDCSSelectFittedSolver<solverType, VC, elementType, degree>::SolverBackendType;
    using Problem =
        typename TDCSSelectFittedSolver<solverType, VC, elementType, degree>::ProblemType;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using ElementSearch = KDTreeElementSearch<typename VC::GridView>;
  };

  template <int dim, ElementType elementType, FittedSolverType solverType, int degree,
            bool geometryAdaption = false>
  class FittedTDCSDriver : public TDCSDriverInterface<dim>
  {
  public:
    using Traits = FittedTDCSDriverTraits<dim, elementType, solverType, degree, geometryAdaption>;

    explicit FittedTDCSDriver(const PatchSet<double, dim>& patchSet,
                              const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : FittedTDCSDriver(FittedDriverData<dim>{}, patchSet, config, dataTree)
    {
    }

    explicit FittedTDCSDriver(const FittedDriverData<dim>& data,
                              const PatchSet<double, dim>& patchSet,
                              const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : config_(config)
        , volumeConductorStorage_(data, config.sub("volume_conductor"),
                                  dataTree.sub("volume_conductor"))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(
              volumeConductorStorage_.get()->gridView()))
        , problem_(
              std::make_shared<typename Traits::Problem>(volumeConductorStorage_.get(), patchSet))
        , solver_(std::make_shared<typename Traits::Solver>(
              volumeConductorStorage_.get(), elementSearch_, problem_,
              config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree()))
        , solverBackend_(std::make_shared<typename Traits::SolverBackend>(
              solver_, config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree()))
    {
    }
    virtual std::unique_ptr<DenseMatrix<double>> CenterEvaluation(const Function& solution)
    {}
    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }

    virtual void solveTDCSForward(Function& solution, const Dune::ParameterTree& config,
                                  DataTree dataTree = DataTree()) override
    {
      solver_->solve(solverBackend_->get(), solution.cast<typename Traits::DomainDOFVector>(),
                     config, dataTree);
    }

    virtual void write(const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        VTKWriter<typename Traits::VC> writer(volumeConductorStorage_.get(),
                                              config.get<unsigned int>("subsampling", degree - 1));
        writer.addCellData(std::make_shared<duneuro::FittedTensorNormFunctor<typename Traits::VC>>(
            volumeConductorStorage_.get()));

#if HAVE_EIGEN
        if (config.get("anisotropy.enable", false)) {
          for (unsigned int i = 0; i < dim; ++i) {
            writer.addCellData(std::make_shared<duneuro::FittedTensorFunctor<typename Traits::VC>>(
                volumeConductorStorage_.get(), i));
          }
        }
#endif

        writer.write(config.get<std::string>("filename"), dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual void write(const Function& function, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        VTKWriter<typename Traits::VC> writer(volumeConductorStorage_.get(),
                                              config.get<unsigned int>("subsampling", degree - 1));
        auto gradient_type = config.get<std::string>("gradient.type", "vertex");
        auto potential_type = config.get<std::string>("potential.type", "vertex");

        if (gradient_type == "vertex") {
          writer.addVertexDataGradient(
              *solver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "gradient_potential");

        } else {
          writer.addCellDataGradient(
              *solver_,
              Dune::stackobject_to_shared_ptr(function.cast<typename Traits::DomainDOFVector>()),
              "gradient_potential");
        }
        if (potential_type == "vertex") {
          writer.addVertexData(*solver_, Dune::stackobject_to_shared_ptr(
                                             function.cast<typename Traits::DomainDOFVector>()),
                               "potential");

        } else {
          writer.addCellData(*solver_, Dune::stackobject_to_shared_ptr(
                                           function.cast<typename Traits::DomainDOFVector>()),
                             "potential");
        }
        writer.addCellData(std::make_shared<duneuro::FittedTensorNormFunctor<typename Traits::VC>>(
            volumeConductorStorage_.get()));

#if HAVE_EIGEN
        if (config.get("anisotropy.enable", false)) {
          for (unsigned int i = 0; i < dim; ++i) {
            writer.addCellData(std::make_shared<duneuro::FittedTensorFunctor<typename Traits::VC>>(
                volumeConductorStorage_.get(), i));
          }
        }
#endif

        writer.write(config.get<std::string>("filename"), dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }
        virtual void
    setElectrodes(const std::vector<typename TDCSDriverInterface<dim>::CoordinateType>& electrodes,
                  const Dune::ParameterTree& config) override
                  {}

  private:
    Dune::ParameterTree config_;
    typename Traits::VCStorage volumeConductorStorage_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Problem> problem_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<typename Traits::SolverBackend> solverBackend_;
  };
}

#endif // DUNEURO_FITTED_TDCS_DRIVER_HH
