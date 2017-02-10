#ifndef DUNEURO_SUBTRACTIONDGRESIDUAL_HH
#define DUNEURO_SUBTRACTIONDGRESIDUAL_HH

#include <dune/common/parametertree.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_default_parameter.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>

namespace duneuro
{
  template <class VC, class FS, class V>
  class SubtractionSourceModel : public SourceModelBase<typename FS::GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename FS::GFS::Traits::GridViewType, V>;
    enum { dim = VC::dim };
    using Problem = SubtractionDGDefaultParameter<typename FS::GFS::Traits::GridViewType,
                                                  typename V::field_type, VC>;
    using EdgeNormProvider = FaceBasedEdgeNormProvider;
    using LOP = SubtractionDG<Problem, EdgeNormProvider>;
    using DOF = typename FS::DOF;
    using AS = Dune::PDELab::GalerkinGlobalAssembler<FS, LOP, Dune::SolverCategory::sequential>;
    using ElementType = typename BaseT::ElementType;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using SearchType = typename BaseT::SearchType;

    SubtractionSourceModel(std::shared_ptr<VC> volumeConductor, const FS& fs,
                           std::shared_ptr<SearchType> search, const Dune::ParameterTree& config)
        : BaseT(search)
        , problem_(volumeConductor->gridView(), volumeConductor)
        , edgeNormProvider_()
        , lop_(problem_, edgeNormProvider_, ConvectionDiffusion_DG_Scheme::SIPG,
               ConvectionDiffusion_DG_Weights::weightsOn, 1.0,
               config.get<unsigned int>("intorderadd"), config.get<unsigned int>("intorderadd_lb"))
        , x_(fs.getGFS(), 0.0)
        , res_(fs.getGFS(), 0.0)
        , interp_(fs.getGFS(), 0.0)
        , assembler_(fs, lop_, 1)
    {
    }

    virtual void bind(const typename BaseT::DipoleType& dipole,
                      DataTree dataTree = DataTree()) override
    {
      BaseT::bind(dipole);
      problem_.bind(this->dipoleElement(), this->localDipolePosition(), this->dipole().moment());
    }

    virtual void assembleRightHandSide(VectorType& vector) const override
    {
      x_ = 0.0;
      assembler_->residual(x_, vector);
      vector *= -1.0;
    }

    virtual void postProcessSolution(VectorType& vector) const override
    {
      interp_ = 0.0;
      Dune::PDELab::interpolate(problem_.get_u_infty(), assembler_->trialGridFunctionSpace(),
                                interp_);
      vector += interp_;
    }

    virtual void
    postProcessSolution(const std::vector<CoordinateType>& electrodes,
                        std::vector<typename VectorType::field_type>& vector) const override
    {
      assert(electrodes.size() == vector.size());
      Dune::FieldVector<typename Problem::Traits::RangeFieldType, 1> result;
      for (unsigned int i = 0; i < electrodes.size(); ++i) {
        problem_.get_u_infty().evaluateGlobal(electrodes[i], result);
        vector[i] += result;
      }
    }

  private:
    Problem problem_;
    EdgeNormProvider edgeNormProvider_;
    LOP lop_;
    mutable DOF x_;
    mutable DOF res_;
    mutable DOF interp_;
    mutable AS assembler_;
  };
}

#endif // DUNEURO_SUBTRACTIONDGRESIDUAL_HH
