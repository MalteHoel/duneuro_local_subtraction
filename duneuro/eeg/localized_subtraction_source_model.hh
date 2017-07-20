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
#include <duneuro/common/subset_entityset.hh>
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

    LocalizedSubtractionSourceModel(std::shared_ptr<VC> volumeConductor,
                                    std::shared_ptr<const FS> fs,
                                    std::shared_ptr<SearchType> search,
                                    const Dune::ParameterTree& config,
                                    const Dune::ParameterTree& solverConfig)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , functionSpace_(fs)
        , elementNeighborhoodMap_(std::make_shared<ElementNeighborhoodMap<typename VC::GridView>>(
              volumeConductor_->gridView()))
        , edgeNormProvider_(solverConfig.get<std::string>("edge_norm_type"), 1.0)
        , config_(config)
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
      lop_ = std::make_shared<LOP>(*problem_, weights_, config_.get<unsigned int>("intorderadd"),
                                   config_.get<unsigned int>("intorderadd_lb"));
      subFS_ = std::make_shared<SUBFS>(subVolumeConductor_);
      dataTree.set("sub_dofs", subFS_->getGFS().size());
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

      SubLFS sublfs(subFS_->getGFS());
      SubLFSCache subcache(sublfs);
      HostLFS hostlfs(functionSpace_->getGFS());
      HostLFSCache hostcache(hostlfs);
      for (const auto& e : Dune::elements(subVolumeConductor_->entitySet())) {
        sublfs.bind(e);
        subcache.update();
        hostlfs.bind(e);
        hostcache.update();
        for (unsigned int i = 0; i < hostcache.size(); ++i) {
          vector[hostcache.containerIndex(i)] += (*x_)[subcache.containerIndex(i)];
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
    std::shared_ptr<const FS> functionSpace_;
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

      SubLFS sublfs(subFS_->getGFS());
      SubLFSCache subcache(sublfs);
      HostLFS hostlfs(functionSpace_->getGFS());
      HostLFSCache hostcache(hostlfs);
      for (const auto& e : Dune::elements(subVolumeConductor_->entitySet())) {
        sublfs.bind(e);
        subcache.update();
        hostlfs.bind(e);
        hostcache.update();
        for (unsigned int i = 0; i < hostcache.size(); ++i) {
          vector[hostcache.containerIndex(i)] = (*r_)[subcache.containerIndex(i)];
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

      SubLFS sublfs(subFS_->getGFS());
      SubLFSCache subcache(sublfs);
      HostLFS hostlfs_inside(functionSpace_->getGFS());
      HostLFSCache hostcache_inside(hostlfs_inside);
      HostLFS hostlfs_outside(functionSpace_->getGFS());
      HostLFSCache hostcache_outside(hostlfs_outside);
      for (const auto& is : patchBoundaryIntersections_) {
        const auto& geo = is.geometry();

        // retrieve and bind inside
        const auto& inside = is.inside();
        hostlfs_inside.bind(inside);
        hostcache_inside.update();
        const auto& A_s = volumeConductor_->tensor(inside);

        // retrieve and bind outside
        const auto& outside = is.outside();
        hostlfs_outside.bind(outside);
        hostcache_outside.update();
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

        const int order_s = FESwitch::basis(hostlfs_inside.finiteElement()).order();
        const int order_n = FESwitch::basis(hostlfs_outside.finiteElement()).order();

        const int degree = std::max(order_s, order_n);
        const RF penalty_factor =
            (penalty_ / h_F) * harmonic_average * degree * (degree + VC::dim - 1);

        const int intorder = intorderadd_lb_ + 2 * degree;

        std::vector<RangeType> phi_s(hostcache_inside.size());
        std::vector<RangeType> phi_n(hostcache_outside.size());
        std::vector<Dune::FieldMatrix<RF, 1, dim>> gradpsi_s(hostcache_inside.size());
        std::vector<Dune::FieldMatrix<RF, 1, dim>> gradpsi_n(hostcache_outside.size());
        std::vector<Dune::FieldVector<RF, dim>> Agradpsi_s(hostcache_inside.size());
        std::vector<Dune::FieldVector<RF, dim>> Agradpsi_n(hostcache_outside.size());
        const auto& rule = Dune::QuadratureRules<DF, VC::dim - 1>::rule(geo.type(), intorder);
        for (const auto& qp : rule) {
          auto qp_inside = is.geometryInInside().global(qp.position());
          auto qp_outside = is.geometryInOutside().global(qp.position());
          // evaluate basis function and their gradients
          FESwitch::basis(hostlfs_inside.finiteElement()).evaluateFunction(qp_inside, phi_s);
          BasisSwitch::gradient(FESwitch::basis(hostlfs_inside.finiteElement()), inside.geometry(),
                                qp_inside, gradpsi_s);
          FESwitch::basis(hostlfs_outside.finiteElement()).evaluateFunction(qp_outside, phi_n);
          BasisSwitch::gradient(FESwitch::basis(hostlfs_outside.finiteElement()),
                                outside.geometry(), qp_outside, gradpsi_n);
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
          for (unsigned int i = 0; i < hostcache_inside.size(); i++) {
            auto index = hostcache_inside.containerIndex(i);
            vector[index] += phi_s[i] * term1;
            vector[index] += -phi_s[i] * omega_n * term1;
            vector[index] += (Agradpsi_s[i] * n_F) * term2 * omega_s; // symmetry term
            vector[index] += -phi_s[i] * term3; // penalty term
          }
          for (unsigned int i = 0; i < hostcache_outside.size(); i++) {
            auto index = hostcache_outside.containerIndex(i);
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