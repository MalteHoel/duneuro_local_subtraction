#ifndef DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH
#define DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH

#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/meg/meg_local_operator.hh>

namespace duneuro
{
  template <class VC, class FluxFS>
  class MEGIntegralAssembler
  {
  public:
    using LOP = MEGLocalOperator<VC, typename FluxFS::FEM>;
    using Assembler =
        Dune::PDELab::GalerkinGlobalAssembler<FluxFS, LOP, Dune::SolverCategory::sequential>;
    using DOF = typename FluxFS::DOF;
    using DomainType = typename LOP::DomainType;

    MEGIntegralAssembler(std::shared_ptr<const VC> vc, std::shared_ptr<const FluxFS> fs,
                         const Dune::ParameterTree& config)
        : fluxFunctionSpace_(fs)
        , dummyx_(fluxFunctionSpace_->getGFS(), 0.0)
        , lop_(vc, config)
        , assembler_(*fluxFunctionSpace_, lop_, 2 * VC::dim + 1)
    {
    }

    void assemble(DOF& dof)
    {
      dof = 0.0;
      assembler_->residual(dummyx_, dof);
    }

    void bind(const DomainType& sensor, const DomainType& projection)
    {
      lop_.bind(sensor, projection);
    }

  private:
    std::shared_ptr<const FluxFS> fluxFunctionSpace_;
    DOF dummyx_;
    LOP lop_;
    Assembler assembler_;
  };

  template <class VC, class FS>
  class CachedIntegralAssembler
  {
  public:
    using Assembler = MEGIntegralAssembler<VC, FS>;
    using DOF = typename FS::DOF;
    using DomainType = Dune::FieldVector<typename VC::ctype, VC::dim>;

    CachedIntegralAssembler(std::shared_ptr<const VC> volumeConductor,
                            std::shared_ptr<const FS> functionSpace,
                            const std::vector<DomainType>& coils,
                            const std::vector<std::vector<DomainType>>& projections,
                            const Dune::ParameterTree& config)
        : cached_(config.get("cache.enable", true))
        , assembler_(volumeConductor, functionSpace, config)
        , coils_(coils)
        , projections_(projections)
        , config_(config)
    {
      int verbose = config.get("verbose", 0);
      if (cached_) {
        allocateCachedDOFs(functionSpace);
        assembleCachedDOFs(assembler_, verbose);
      }
    }

    void assemble(unsigned int coilIndex, unsigned int projectionIndex, DOF& dof)
    {
      if (cached_) {
        dof = *cachedDofs_[coilIndex][projectionIndex];
      } else {
#if HAVE_TBB
        assemble(assembler_.local(), coilIndex, projectionIndex, dof);
#else
        assemble(assembler_, coilIndex, projectionIndex, dof);
#endif
      }
    }

  private:
    bool cached_;
#if HAVE_TBB
    tbb::enumerable_thread_specific<Assembler> assembler_;
#else
    Assembler assembler_;
#endif
    std::vector<std::vector<std::shared_ptr<DOF>>> cachedDofs_;
    std::vector<DomainType> coils_;
    std::vector<std::vector<DomainType>> projections_;
    Dune::ParameterTree config_;

    void assemble(Assembler& assembler, unsigned int coilIndex, unsigned int projectionIndex,
                  DOF& dof) const
    {
      assembler.bind(coils_[coilIndex], projections_[coilIndex][projectionIndex]);
      assembler.assemble(dof);
    }

    void allocateCachedDOFs(std::shared_ptr<const FS> functionSpace)
    {
      cachedDofs_.resize(coils_.size());
      for (unsigned int coil = 0; coil < coils_.size(); ++coil) {
        std::vector<std::shared_ptr<DOF>> dofs_;
        for (unsigned int projection = 0; projection < projections_[coil].size(); ++projection) {
          dofs_.push_back(std::make_shared<DOF>(functionSpace->getGFS(), 0.0));
        }
        cachedDofs_[coil] = dofs_;
      }
    }

    void assembleCachedDOFs(Assembler& assembler, int verbose)
    {
      for (unsigned int coil = 0; coil < coils_.size(); ++coil) {
        for (unsigned int projection = 0; projection < projections_[coil].size(); ++projection) {
          Dune::Timer timer;
          assemble(assembler, coil, projection, *cachedDofs_[coil][projection]);
          if (verbose > 0) {
            std::cout << "assembled integral for coil " << coil << " projection " << projection
                      << ": " << timer.elapsed() << " s" << std::endl;
          }
        }
      }
    }

#if HAVE_TBB
    void assembleCachedDOFs(tbb::enumerable_thread_specific<Assembler>& assembler, int verbose)
    {
      tbb::task_scheduler_init init(config_.hasKey("numberOfThreads") ?
                                        config_.get<std::size_t>("numberOfThreads") :
                                        tbb::task_scheduler_init::automatic);
      auto grainSize = config_.get<int>("grainSize", 16);
      // split coils into blocks of at most grainSize entries and assemble in parallel
      tbb::parallel_for(
          tbb::blocked_range<std::size_t>(0, coils_.size(), grainSize),
          [&](const tbb::blocked_range<std::size_t>& range) {
            for (unsigned int coil = range.begin(); coil < range.end(); ++coil) {
              for (unsigned int projection = 0; projection < projections_[coil].size();
                   ++projection) {
                Dune::Timer timer;
                assemble(assembler.local(), coil, projection, *cachedDofs_[coil][projection]);
                if (verbose > 0) {
                  std::cout << "assembled integral for coil " << coil << " projection "
                            << projection << ": " << timer.elapsed() << " s" << std::endl;
                }
              }
            }
          });
    }
#endif
  };
}

#endif // DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH
