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
    {
      int verbose = config.get("verbose", 0);
      if (cached_) {
        for (unsigned int coil = 0; coil < coils.size(); ++coil) {
          std::vector<std::shared_ptr<DOF>> dofs_;
          for (unsigned int projection = 0; projection < projections[coil].size(); ++projection) {
            Dune::Timer timer;
            auto dof = std::make_shared<DOF>(functionSpace->getGFS(), 0.0);
            assembler_.bind(coils[coil], projections[coil][projection]);
            assembler_.assemble(*dof);
            dofs_.push_back(dof);
            if (verbose > 0) {
              std::cout << "assembled integral for coil " << coil << " projection " << projection
                        << ": " << timer.elapsed() << " s" << std::endl;
            }
          }
          cachedDofs_.push_back(dofs_);
        }
      } else {
        coils_ = coils;
        projections_ = projections;
        uncachedDof_ = std::make_shared<DOF>(functionSpace->getGFS(), 0.0);
      }
    }

    const DOF& assemble(unsigned int coilIndex, unsigned int projectionIndex)
    {
      if (cached_) {
        return *cachedDofs_[coilIndex][projectionIndex];
      } else {
        assembler_.bind(coils_[coilIndex], projections_[coilIndex][projectionIndex]);
        assembler_.assemble(*uncachedDof_);
        return *uncachedDof_;
      }
    }

  private:
    bool cached_;
    std::shared_ptr<DOF> uncachedDof_;
    Assembler assembler_;
    std::vector<std::vector<std::shared_ptr<DOF>>> cachedDofs_;
    std::vector<DomainType> coils_;
    std::vector<std::vector<DomainType>> projections_;
  };
}

#endif // DUNEURO_MEG_INTEGRAL_ASSEMBLER_HH
