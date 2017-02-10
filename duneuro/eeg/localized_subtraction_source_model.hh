#ifndef DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH
#define DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH

#include <dune/common/parametertree.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>

#include <dune/subgrid/subgrid.hh>

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/common/logged_timer.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_default_parameter.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>
#include <duneuro/io/vtk_writer.hh>

namespace duneuro
{
  template <class SubGrid, class HostFS>
  class SubFunctionSpace
  {
  public:
    using Grid = SubGrid;
    using GV = typename Grid::LeafGridView;
    using FEM = typename HostFS::FEM;
    using CONB = typename HostFS::CONB;
    using CON = typename HostFS::CON;
    using VBE = typename HostFS::VBE;
    using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
    using DOF = Dune::PDELab::Backend::Vector<GFS, typename HostFS::NT>;
    using DGF = Dune::PDELab::DiscreteGridFunction<GFS, DOF>;
    using CC = typename GFS::template ConstraintsContainer<typename HostFS::NT>::Type;
    using NT = typename HostFS::NT;

    explicit SubFunctionSpace(const GV& gridView)
        : gv(gridView)
        , femp(std::make_shared<FEM>())
        , gfsp(std::make_shared<GFS>(gv, *femp))
        , ccp(std::make_shared<CC>())
    {
      gfsp->update();
    }

    FEM& getFEM()
    {
      return *femp;
    }
    const FEM& getFEM() const
    {
      return *femp;
    }

    // return gfs reference
    GFS& getGFS()
    {
      return *gfsp;
    }

    // return gfs reference const version
    const GFS& getGFS() const
    {
      return *gfsp;
    }
    // return gfs reference
    CC& getCC()
    {
      return *ccp;
    }

    // return gfs reference const version
    const CC& getCC() const
    {
      return *ccp;
    }

    template <class BCTYPE>
    void assembleConstraints(const BCTYPE& bctype)
    {
      ccp->clear();
      constraints(bctype, *gfsp, *ccp);
    }

    void clearConstraints()
    {
      ccp->clear();
    }

    void setConstrainedDOFS(DOF& x, typename HostFS::NT nt) const
    {
      set_constrained_dofs(*ccp, nt, x);
      conb.make_consistent(*gfsp, x);
    }

    void setNonConstrainedDOFS(DOF& x, typename HostFS::NT nt) const
    {
      set_nonconstrained_dofs(*ccp, nt, x);
      conb.make_consistent(*gfsp, x);
    }

    void copyConstrainedDOFS(const DOF& xin, DOF& xout) const
    {
      copy_constrained_dofs(*ccp, xin, xout);
      conb.make_consistent(*gfsp, xout);
    }

    void copyNonConstrainedDOFS(const DOF& xin, DOF& xout) const
    {
      copy_nonconstrained_dofs(*ccp, xin, xout);
      conb.make_consistent(*gfsp, xout);
    }

  private:
    GV gv; // need this object here because FEM and GFS store a const reference !!
    CONB conb;
    std::shared_ptr<FEM> femp;
    std::shared_ptr<GFS> gfsp;
    std::shared_ptr<CC> ccp;
  };

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
    using SubGrid = Dune::SubGrid<dim, typename VC::GridType>;
    using HostGridView = typename VC::GridView;
    using SubGridView = typename SubGrid::LeafGridView;
    using SubVolumeConductor = VolumeConductor<SubGrid>;
    using Problem =
        SubtractionDGDefaultParameter<SubGridView, typename V::field_type, SubVolumeConductor>;
    using EdgeNormProvider = MultiEdgeNormProvider;
    using LOP = SubtractionDG<Problem, EdgeNormProvider>;
    using SUBFS = SubFunctionSpace<SubGrid, FS>;
    using DOF = typename SUBFS::DOF;
    using AS = Dune::PDELab::GalerkinGlobalAssembler<SUBFS, LOP, Dune::SolverCategory::sequential>;
    using SubLFS = Dune::PDELab::LocalFunctionSpace<typename SUBFS::GFS>;
    using SubLFSCache = Dune::PDELab::LFSIndexCache<SubLFS>;
    using HostLFS = Dune::PDELab::LocalFunctionSpace<typename FS::GFS>;
    using HostLFSCache = Dune::PDELab::LFSIndexCache<HostLFS>;
    using HostProblem = SubtractionDGDefaultParameter<HostGridView, typename V::field_type, VC>;

    LocalizedSubtractionSourceModel(std::shared_ptr<VC> volumeConductor, const FS& fs,
                                    std::shared_ptr<SearchType> search,
                                    const Dune::ParameterTree& config)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , elementNeighborhoodMap_(std::make_shared<ElementNeighborhoodMap<typename VC::GridView>>(
              volumeConductor_->gridView()))
        , edgeNormProvider_(config.get<std::string>("edge_norm_type"), 1.0)
        , config_(config)
        , lfsInside_(fs.getGFS())
        , lfsCacheInside_(lfsInside_)
        , lfsOutside_(fs.getGFS())
        , lfsCacheOutside_(lfsOutside_)
        , intorderadd_lb_(config.get<unsigned int>("intorderadd_lb"))
        , weights_(config.get<bool>("weights") ? ConvectionDiffusion_DG_Weights::weightsOn :
                                                 ConvectionDiffusion_DG_Weights::weightsOff)
        , penalty_(config.get<double>("penalty"))
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

      // create subgrid for selected elements
      std::unique_ptr<SubGrid> subGrid;
      if (!subVolumeConductor_) {
        subGrid = Dune::Std::make_unique<SubGrid>(volumeConductor_->grid());
      } else {
        subGrid.reset(subVolumeConductor_->releaseGrid());
      }
      timer.lap("sub_grid_init");
      subGrid->createBegin();
      timer.lap("sub_grid_create_begin");
      for (const auto& element : elementPatch->elements()) {
        subGrid->insert(element);
      }
      timer.lap("sub_grid_insert_elements");
      subGrid->createEnd(false);
      timer.lap("sub_grid_create_end");
      auto subGridView = subGrid->leafGridView();
      timer.lap("sub_grid_create_view");

      // extract conductivity tensors to create a local volume conductor
      Dune::SingleCodimSingleGeomTypeMapper<SubGridView, 0> subGridElementMapper(subGridView);
      std::vector<typename VC::TensorType> tensors(subGridElementMapper.size());
      for (const auto& subElement : Dune::elements(subGridView)) {
        tensors[subGridElementMapper.index(subElement)] =
            volumeConductor_->tensor(subGrid->template getHostEntity<0>(subElement));
      }
      timer.lap("extract_sub_tensors");

      // create sub grid volume conductor
      subVolumeConductor_ = std::make_shared<SubVolumeConductor>(
          std::move(subGrid),
          Dune::Std::make_unique<typename SubVolumeConductor::MappingType>(
              DirectEntityMapping<SubGridView, typename VC::TensorType>(subGridView, tensors)));
      timer.lap("sub_volume_conductor");

      if (config_.hasSub("vtk")) {
        if (config_.get("vtk.enable", false)) {
          VTKWriter<SubVolumeConductor, 1> writer(subVolumeConductor_,
                                                  config_.get("vtk.subsampling", 0));
          writer.write(config_.get<std::string>("vtk.filename"));
          timer.lap("vtk");
        }
      }

      hostProblem_ = std::make_shared<HostProblem>(volumeConductor_->gridView(), volumeConductor_);
      hostProblem_->bind(this->dipoleElement(), this->localDipolePosition(),
                         this->dipole().moment());
      problem_ = std::make_shared<Problem>(subVolumeConductor_->gridView(), subVolumeConductor_);
      problem_->bind(
          subVolumeConductor_->grid().template getSubGridEntity<0>(this->dipoleElement()),
          this->localDipolePosition(), this->dipole().moment());
      lop_ = std::make_shared<LOP>(
          *problem_, edgeNormProvider_, ConvectionDiffusion_DG_Scheme::SIPG,
          ConvectionDiffusion_DG_Weights::weightsOn, 1.0, config_.get<unsigned int>("intorderadd"),
          config_.get<unsigned int>("intorderadd_lb"));
      subFS_ = std::make_shared<SUBFS>(subVolumeConductor_->gridView());
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

      Dune::SingleCodimSingleGeomTypeMapper<SubGridView, 0> subGridElementMapper(
          subVolumeConductor_->gridView());
      Dune::SingleCodimSingleGeomTypeMapper<HostGridView, 0> hostGridElementMapper(
          volumeConductor_->gridView());
      using Dune::PDELab::Backend::native;
      for (const auto& e : Dune::elements(subVolumeConductor_->gridView())) {
        sublfsInside_->bind(e);
        sublfsCacheInside_->update();
        const auto& hostEntity = subVolumeConductor_->grid().template getHostEntity<0>(e);
        lfsInside_.bind(hostEntity);
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
    ConvectionDiffusion_DG_Weights::Type weights_;
    double penalty_;

    void assembleLocalDefaultSubtraction(VectorType& vector) const
    {
      *x_ = 0.0;
      *r_ = 0.0;
      (*assembler_)->residual(*x_, *r_);
      *r_ *= -1.;

      Dune::SingleCodimSingleGeomTypeMapper<SubGridView, 0> subGridElementMapper(
          subVolumeConductor_->gridView());
      Dune::SingleCodimSingleGeomTypeMapper<HostGridView, 0> hostGridElementMapper(
          volumeConductor_->gridView());
      for (const auto& e : Dune::elements(subVolumeConductor_->gridView())) {
        sublfsInside_->bind(e);
        sublfsCacheInside_->update();
        const auto& hostEntity = subVolumeConductor_->grid().template getHostEntity<0>(e);
        lfsInside_.bind(hostEntity);
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

        const auto& rule = Dune::QuadratureRules<DF, VC::dim - 1>::rule(geo.type(), intorder);
        for (const auto& qp : rule) {
          /** evaluate test and ansatz functions **/
          std::vector<RangeType> phi_s(lfsCacheInside_.size());
          FESwitch::basis(lfsInside_.finiteElement())
              .evaluateFunction(is.geometryInInside().global(qp.position()), phi_s);
          std::vector<RangeType> phi_n(lfsCacheOutside_.size());
          FESwitch::basis(lfsOutside_.finiteElement())
              .evaluateFunction(is.geometryInOutside().global(qp.position()), phi_n);

          RF factor = qp.weight() * geo.integrationElement(qp.position());

          auto global = geo.global(qp.position());

          Dune::FieldVector<RF, VC::dim> A_s_graduinfty;
          auto graduinfty = hostProblem_->get_grad_u_infty(global);
          A_s.mv(graduinfty, A_s_graduinfty);
          auto n_F_A_s_graduinfty = n_F * A_s_graduinfty;
          auto term1 = n_F_A_s_graduinfty * factor;
          auto cf_flux = config_.get<double>("cf_flux");
          for (unsigned int i = 0; i < lfsCacheInside_.size(); i++) {
            vector[lfsCacheInside_.containerIndex(i)] += phi_s[i] * term1;
            vector[lfsCacheInside_.containerIndex(i)] += -phi_s[i] * omega_n * term1 * cf_flux;
          }
          for (unsigned int i = 0; i < lfsCacheOutside_.size(); i++) {
            vector[lfsCacheOutside_.containerIndex(i)] += -phi_n[i] * omega_s * term1 * cf_flux;
          }
          auto cf_jump = config_.get<double>("cf_jump");
          auto uinfty = hostProblem_->get_u_infty(global);
          auto term2 = factor * penalty_factor * uinfty;
          for (unsigned int i = 0; i < lfsCacheInside_.size(); i++) {
            vector[lfsCacheInside_.containerIndex(i)] += -phi_s[i] * term2 * cf_jump;
          }
          for (unsigned int i = 0; i < lfsCacheOutside_.size(); i++) {
            vector[lfsCacheOutside_.containerIndex(i)] += phi_n[i] * term2 * cf_jump;
          }
        }
      }
    }
  };
}

#endif // DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH
