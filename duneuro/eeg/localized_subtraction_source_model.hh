#ifndef DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH
#define DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH

#include <dune/common/parametertree.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/common/entityset_volume_conductor.hh>
#include <duneuro/common/logged_timer.hh>
#include <duneuro/common/sub_function_space.hh>
#include <duneuro/common/subsetentityset.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_default_parameter.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>

namespace duneuro
{
  template <class VC, class FS, class V>
  class LocalizedSubtractionSourceModel
      : public SourceModelBase<typename FS::GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename FS::GFS::Traits::GridViewType, V>;
    enum { dim = VC::dim };
    using ElementType = typename BaseT::ElementType;
    using CoordinateType = typename BaseT::CoordinateType;
    using VectorType = typename BaseT::VectorType;
    using SearchType = typename BaseT::SearchType;
    using HostGridView = typename VC::GridView;
    using SubEntitySet = SubSetEntitySet<HostGridView>;
    using SubVolumeConductor = EntitySetVolumeConductor<SubEntitySet>;
    using SUBFS = SubFunctionSpace<FS, SubVolumeConductor>;
    using Problem =
        SubtractionDGDefaultParameter<SubEntitySet, typename V::field_type, SubVolumeConductor>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using LOP = SubtractionDG<Problem, EdgeNormProvider>;
    using DOF = typename SUBFS::DOF;
    using AS = Dune::PDELab::GalerkinGlobalAssembler<SUBFS, LOP, Dune::SolverCategory::sequential>;
    using SubLFS = Dune::PDELab::LocalFunctionSpace<typename SUBFS::GFS>;
    using SubLFSCache = Dune::PDELab::LFSIndexCache<SubLFS>;
    using HostLFS = Dune::PDELab::LocalFunctionSpace<typename FS::GFS>;
    using HostLFSCache = Dune::PDELab::LFSIndexCache<HostLFS>;
    using HostProblem = SubtractionDGDefaultParameter<HostGridView, typename V::field_type, VC>;

    LocalizedSubtractionSourceModel(std::shared_ptr<VC> volumeConductor, const FS& fs,
                                    std::shared_ptr<SearchType> search,
                                    const Dune::ParameterTree& config,
                                    const Dune::ParameterTree& solverConfig)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , elementNeighborhoodMap_(std::make_shared<ElementNeighborhoodMap<typename VC::GridView>>(
              volumeConductor_->gridView()))
        , edgeNormProvider_(solverConfig.get<std::string>("edge_norm_type"), 1.0)
        , config_(config)
        , lfsInside_(fs.getGFS())
        , lfsCacheInside_(lfsInside_)
        , lfsOutside_(fs.getGFS())
        , lfsCacheOutside_(lfsOutside_)
        , intorderadd_lb_(config.get<unsigned int>("intorderadd_lb"))
        , scheme_(
              ConvectionDiffusion_DG_Scheme::fromString(solverConfig.get<std::string>("scheme")))
        , weights_(solverConfig.get<bool>("weights") ? ConvectionDiffusion_DG_Weights::weightsOn :
                                                       ConvectionDiffusion_DG_Weights::weightsOff)
        , penalty_(solverConfig.get<double>("penalty"))
    {
    }

    virtual void bind(const typename BaseT::DipoleType& dipole,
                      DataTree dataTree = DataTree()) override
    {
      LoggedTimer timer(dataTree);
      BaseT::bind(dipole, dataTree);
      timer.lap("bind_base");

      // select elements for the local subtraction space
      auto elementPatch = make_element_patch(volumeConductor_, elementNeighborhoodMap_,
                                             this->elementSearch(), dipole.position(), config_);
      timer.lap("make_element_patch");
      dataTree.set("elements", elementPatch->elements().size());

      // extract patch boundary intersection
      patchBoundaryIntersections_ = elementPatch->extractBoundaryIntersections();
      timer.lap("extract_boundary_intersection");

      SubEntitySet subEntitySet(volumeConductor_->gridView(), elementPatch->elements());
      timer.lap("create_sub_entity_set");

      // extract conductivity tensors to create a local volume conductor
      Dune::MultipleCodimMultipleGeomTypeMapper<SubEntitySet, Dune::MCMGElementLayout> mapper(
          subEntitySet);
      std::vector<typename VC::TensorType> tensors(mapper.size());
      for (const auto& subElement : elementPatch->elements()) {
        tensors[mapper.index(subElement)] = volumeConductor_->tensor(subElement);
      }
      timer.lap("extract_sub_tensors");

      // create sub grid volume conductor
      subVolumeConductor_ = std::make_shared<SubVolumeConductor>(subEntitySet, tensors);
      timer.lap("sub_volume_conductor");

      hostProblem_ = std::make_shared<HostProblem>(volumeConductor_->gridView(), volumeConductor_);
      hostProblem_->bind(this->dipoleElement(), this->localDipolePosition(),
                         this->dipole().moment());
      problem_ = std::make_shared<Problem>(subVolumeConductor_->entitySet(), subVolumeConductor_);
      problem_->bind(this->dipoleElement(), this->localDipolePosition(), this->dipole().moment());
      lop_ = std::make_shared<LOP>(*problem_, edgeNormProvider_, scheme_, weights_, penalty_,
                                   config_.get<unsigned int>("intorderadd"),
                                   config_.get<unsigned int>("intorderadd_lb"));
      subFS_ = std::make_shared<SUBFS>(subVolumeConductor_);
      dataTree.set("sub_dofs", subFS_->getGFS().size());
      sublfsInside_ = std::make_shared<SubLFS>(subFS_->getGFS());
      sublfsCacheInside_ = std::make_shared<SubLFSCache>(*sublfsInside_);
      x_ = std::make_shared<DOF>(subFS_->getGFS(), 0.0);
      r_ = std::make_shared<DOF>(subFS_->getGFS(), 0.0);
      assembler_ = std::make_shared<AS>(*subFS_, *lop_, 1);
      timer.lap("sub_problem");
      timer.stop("bind_accumulated");
      // note: maybe invert normal in boundary condition of subtraction operator??
    }

    virtual void assembleRightHandSide(VectorType& vector) const override
    {
      assembleLocalDefaultSubtraction(vector);
      assemblePatchBoundaryTerm(vector);
    }

    virtual void postProcessSolution(VectorType& vector) const override
    {
      *x_ = 0.0;
      Dune::PDELab::interpolate(problem_->get_u_infty(), (*assembler_)->trialGridFunctionSpace(),
                                *x_);

      for (const auto& e : Dune::elements(subVolumeConductor_->entitySet())) {
        sublfsInside_->bind(e);
        sublfsCacheInside_->update();
        lfsInside_.bind(e);
        lfsCacheInside_.update();
        for (unsigned int i = 0; i < lfsInside_.size(); ++i) {
          vector[lfsCacheInside_.containerIndex(i)] += (*x_)[sublfsCacheInside_->containerIndex(i)];
        }
      }
    }

    virtual void
    postProcessSolution(const std::vector<CoordinateType>& electrodes,
                        std::vector<typename VectorType::field_type>& vector) const override
    {
      // note: need to check if electrode is within the patch
      // currently assume that the patch does not touch the boundary
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    std::shared_ptr<ElementNeighborhoodMap<typename VC::GridView>> elementNeighborhoodMap_;
    std::shared_ptr<SubVolumeConductor> subVolumeConductor_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<HostProblem> hostProblem_;
    EdgeNormProvider edgeNormProvider_;
    std::shared_ptr<LOP> lop_;
    std::shared_ptr<SUBFS> subFS_;
    std::shared_ptr<AS> assembler_;
    std::shared_ptr<DOF> x_;
    std::shared_ptr<DOF> r_;
    Dune::ParameterTree config_;
    std::vector<typename HostGridView::Intersection> patchBoundaryIntersections_;
    mutable HostLFS lfsInside_;
    mutable HostLFSCache lfsCacheInside_;
    mutable HostLFS lfsOutside_;
    mutable HostLFSCache lfsCacheOutside_;
    mutable std::shared_ptr<SubLFS> sublfsInside_;
    mutable std::shared_ptr<SubLFSCache> sublfsCacheInside_;
    unsigned int intorderadd_lb_;
    ConvectionDiffusion_DG_Scheme::Type scheme_;
    ConvectionDiffusion_DG_Weights::Type weights_;
    double penalty_;

    void assembleLocalDefaultSubtraction(VectorType& vector) const
    {
      *x_ = 0.0;
      *r_ = 0.0;
      (*assembler_)->residual(*x_, *r_);
      *r_ *= -1.;

      for (const auto& e : Dune::elements(subVolumeConductor_->entitySet())) {
        sublfsInside_->bind(e);
        sublfsCacheInside_->update();
        lfsInside_.bind(e);
        lfsCacheInside_.update();
        for (unsigned int i = 0; i < lfsInside_.size(); ++i) {
          vector[lfsCacheInside_.containerIndex(i)] = (*r_)[sublfsCacheInside_->containerIndex(i)];
        }
      }
    }

    void assemblePatchBoundaryTerm(VectorType& vector) const
    {
      using FESwitch =
          Dune::FiniteElementInterfaceSwitch<typename HostLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using DF = typename BasisSwitch::DomainField;
      using RF = typename BasisSwitch::RangeField;
      using RangeType = typename BasisSwitch::Range;

      for (const auto& is : patchBoundaryIntersections_) {
        const auto& geo = is.geometry();

        // retrieve and bind inside
        const auto& inside = is.inside();
        lfsInside_.bind(inside);
        lfsCacheInside_.update();
        const auto& A_s = volumeConductor_->tensor(inside);

        // retrieve and bind outside
        const auto& outside = is.outside();
        lfsOutside_.bind(outside);
        lfsCacheOutside_.update();
        const auto& A_n = volumeConductor_->tensor(outside);

        const auto& n_F = is.centerUnitOuterNormal();
        RF omega_s;
        RF omega_n;
        RF harmonic_average;
        if (weights_ == ConvectionDiffusion_DG_Weights::weightsOn) {
          Dune::FieldVector<RF, dim> An_F_s;
          A_s.mv(n_F, An_F_s);
          Dune::FieldVector<RF, dim> An_F_n;
          A_n.mv(n_F, An_F_n);
          RF delta_s = (An_F_s * n_F);
          RF delta_n = (An_F_n * n_F);
          omega_s = delta_n / (delta_s + delta_n + 1e-20);
          omega_n = delta_s / (delta_s + delta_n + 1e-20);
          harmonic_average = 2.0 * delta_s * delta_n / (delta_s + delta_n + 1e-20);
        } else {
          omega_s = omega_n = 0.5;
          harmonic_average = 1.0;
        }

        // note: edgenorm provider needs the intersectiongeometry interface
        RF h_F;
        edgeNormProvider_.edgeNorm(
            Dune::PDELab::IntersectionGeometry<typename HostGridView::Intersection>(is, 0), h_F);

        const int order_s = FESwitch::basis(lfsInside_.finiteElement()).order();
        const int order_n = FESwitch::basis(lfsOutside_.finiteElement()).order();

        const int degree = std::max(order_s, order_n);
        const RF penalty_factor =
            (penalty_ / h_F) * harmonic_average * degree * (degree + VC::dim - 1);

        const int intorder = intorderadd_lb_ + 2 * degree;

        std::vector<RangeType> phi_s(lfsCacheInside_.size());
        std::vector<RangeType> phi_n(lfsCacheOutside_.size());
        std::vector<Dune::FieldMatrix<RF, 1, dim>> gradpsi_s(lfsCacheInside_.size());
        std::vector<Dune::FieldMatrix<RF, 1, dim>> gradpsi_n(lfsCacheOutside_.size());
        std::vector<Dune::FieldVector<RF, dim>> Agradpsi_s(lfsCacheInside_.size());
        std::vector<Dune::FieldVector<RF, dim>> Agradpsi_n(lfsCacheOutside_.size());
        const auto& rule = Dune::QuadratureRules<DF, VC::dim - 1>::rule(geo.type(), intorder);
        for (const auto& qp : rule) {
          auto qp_inside = is.geometryInInside().global(qp.position());
          auto qp_outside = is.geometryInOutside().global(qp.position());
          // evaluate basis function and their gradients
          FESwitch::basis(lfsInside_.finiteElement()).evaluateFunction(qp_inside, phi_s);
          BasisSwitch::gradient(FESwitch::basis(lfsInside_.finiteElement()), inside.geometry(),
                                qp_inside, gradpsi_s);
          FESwitch::basis(lfsOutside_.finiteElement()).evaluateFunction(qp_outside, phi_n);
          BasisSwitch::gradient(FESwitch::basis(lfsOutside_.finiteElement()), outside.geometry(),
                                qp_outside, gradpsi_n);
          // compute sigma*gradient_psi
          for (unsigned int i = 0; i < gradpsi_s.size(); ++i)
            A_s.mv(gradpsi_s[i][0], Agradpsi_s[i]);
          for (unsigned int i = 0; i < gradpsi_n.size(); ++i)
            A_n.mv(gradpsi_n[i][0], Agradpsi_n[i]);

          RF factor = qp.weight() * geo.integrationElement(qp.position());

          // compute infinity potential and its gradient
          auto global = geo.global(qp.position());
          auto uinfty = hostProblem_->get_u_infty(global);
          auto graduinfty = hostProblem_->get_grad_u_infty(global);
          Dune::FieldVector<RF, VC::dim> A_s_graduinfty;
          A_s.mv(graduinfty, A_s_graduinfty);

          // assemble the integrals
          auto term1 = factor * (n_F * A_s_graduinfty);
          auto term2 = factor * uinfty;
          auto term3 = term2 * penalty_factor;
          for (unsigned int i = 0; i < lfsCacheInside_.size(); i++) {
            auto index = lfsCacheInside_.containerIndex(i);
            vector[index] += phi_s[i] * term1;
            vector[index] += -phi_s[i] * omega_n * term1;
            vector[index] += (Agradpsi_s[i] * n_F) * term2 * omega_s; // symmetry term
            vector[index] += -phi_s[i] * term3; // penalty term
          }
          for (unsigned int i = 0; i < lfsCacheOutside_.size(); i++) {
            auto index = lfsCacheOutside_.containerIndex(i);
            vector[index] += -phi_n[i] * omega_s * term1;
            vector[index] += (Agradpsi_n[i] * n_F) * term2 * omega_n; // symmetry term
            vector[index] += phi_n[i] * term3; // penalty term
          }
        }
      }
    }
  };
}

#endif // DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH
