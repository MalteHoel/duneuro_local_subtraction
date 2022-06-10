#ifndef DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH
#define DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH

#include <type_traits>
#include <unordered_set>

#include <dune/common/parametertree.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>             // include for makeGridViewFunction

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>             // include for the DiscreteGridViewFunction class
#include <dune/pdelab/common/crossproduct.hh>				                    // include for cross product

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/common/edge_norm_provider.hh>
#include <duneuro/common/element_patch.hh>
#include <duneuro/common/entityset_volume_conductor.hh>
#include <duneuro/common/element_patch_assembler.hh>
#include <duneuro/common/logged_timer.hh>
#include <duneuro/common/penalty_flux_weighting.hh>
#include <duneuro/common/sub_function_space.hh>
#include <duneuro/common/subset_entityset.hh>
#include <duneuro/eeg/source_model_interface.hh>
#include <duneuro/eeg/subtraction_dg_default_parameter.hh>
#include <duneuro/eeg/subtraction_dg_operator.hh>
#include <duneuro/eeg/localized_subtraction_dg_local_operator.hh>
#include <duneuro/eeg/localized_subtraction_cg_local_operator.hh>
#include <duneuro/eeg/analytic_utilities.hh>
#include <duneuro/eeg/localized_subtraction_cg_p1_local_operator.hh>
#include <duneuro/common/flags.hh>

namespace duneuro
{
  template <class VC, class FS, class V, ContinuityType continuityType>
  class LocalizedSubtractionSourceModel
      : public SourceModelBase<typename FS::GFS::Traits::GridViewType, V>
  {
  public:
    using BaseT = SourceModelBase<typename FS::GFS::Traits::GridViewType, V>;
    enum { dim = VC::dim };
    using ElementType = typename BaseT::ElementType;
    using CoordinateField = typename ElementType::Geometry::ctype;
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
    using PenaltyFluxWeighting = FittedDynamicPenaltyFluxWeights;
    using LOP = SubtractionDG<FS, Problem, EdgeNormProvider, PenaltyFluxWeighting>;
    using DOF = typename SUBFS::DOF;
    using AS = Dune::PDELab::GalerkinGlobalAssembler<SUBFS, LOP, Dune::SolverCategory::sequential>;
    using SubLFS = Dune::PDELab::LocalFunctionSpace<typename SUBFS::GFS>;
    using SubLFSCache = Dune::PDELab::LFSIndexCache<SubLFS>;
    using HostLFS = Dune::PDELab::LocalFunctionSpace<typename FS::GFS>;
    using HostLFSCache = Dune::PDELab::LFSIndexCache<HostLFS>;
    using HostProblem = SubtractionDGDefaultParameter<HostGridView, typename V::field_type, VC>;
    using DOFVector = Dune::PDELab::Backend::Vector<typename FS::GFS, typename FS::NT>;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename FS::GFS, DOFVector>;
    using LocalFunction = typename DiscreteGridFunction::LocalFunction;
    enum {diffOrder = 1};
    using DerivativeGridFunction = typename Dune::PDELab::DiscreteGridViewFunction<typename DiscreteGridFunction::GridFunctionSpace, DOFVector, diffOrder>;
    using LocalDerivativeFunction = typename DerivativeGridFunction::LocalFunction;
    using Tensor = typename HostProblem::Traits::PermTensorType;

    LocalizedSubtractionSourceModel(std::shared_ptr<const VC> volumeConductor,
                                    std::shared_ptr<const FS> fs,
                                    std::shared_ptr<const SearchType> search,
                                    const Dune::ParameterTree& config,
                                    const Dune::ParameterTree& solverConfig)
        : BaseT(search)
        , volumeConductor_(volumeConductor)
        , functionSpace_(fs)
        , config_(config)
        , patchAssembler_(volumeConductor_,functionSpace_,search,config_)
        // parameters for LOP
        , edgeNormProvider_(solverConfig.get<std::string>("edge_norm_type"), 1.0)
        , weighting_(solverConfig.get<std::string>("weights"))
        , intorderadd_eeg_patch_(config.get<unsigned int>("intorderadd_eeg_patch"))
        , intorderadd_eeg_boundary_(config.get<unsigned int>("intorderadd_eeg_boundary"))
        , intorderadd_eeg_transition_(config.get<unsigned int>("intorderadd_eeg_transition"))
        , intorder_meg_patch_(config.get<unsigned int>("intorder_meg_patch", 0))
        , intorder_meg_boundary_(config.get<unsigned int>("intorder_meg_boundary", 6))
        , intorder_meg_transition_(config.get<unsigned int>("intorder_meg_transition", 5))
        , penalty_(solverConfig.get<double>("penalty"))
        , chiFunctionPtr_(nullptr)
    {
    }

    virtual void bind(const typename BaseT::DipoleType& dipole,
                      DataTree dataTree = DataTree()) override
    {
      LoggedTimer timer(dataTree);
      BaseT::bind(dipole, dataTree);
      timer.lap("bind_base");

      // build patch
      patchAssembler_.bind(dipole.position(), dataTree);
      timer.lap("bind_patch_assembler");
    
      // setup of problem parameters (everything related to the evaluation of u_infinity, its gradient and sigma_infinity)
      hostProblem_ = std::make_shared<HostProblem>(volumeConductor_->gridView(), volumeConductor_);
      hostProblem_->bind(this->dipoleElement(), this->localDipolePosition(),
                         this->dipole().moment());

      if constexpr(continuityType == ContinuityType::discontinuous)
      {
        SubEntitySet subEntitySet(volumeConductor_->gridView(), patchAssembler_.patchElements());
        timer.lap("create_sub_entity_set");

        // extract conductivity tensors to create a local volume conductor
        Dune::MultipleCodimMultipleGeomTypeMapper<SubEntitySet>
            mapper(subEntitySet, Dune::mcmgElementLayout());
        std::vector<typename VC::TensorType> tensors(mapper.size());
        for (const auto& subElement : patchAssembler_.patchElements()) {
          tensors[mapper.index(subElement)] = volumeConductor_->tensor(subElement);
        }
        timer.lap("extract_sub_tensors");

        // create sub grid volume conductor
        subVolumeConductor_ = std::make_shared<SubVolumeConductor>(subEntitySet, tensors);
        timer.lap("sub_volume_conductor");

        problem_ = std::make_shared<Problem>(subVolumeConductor_->entitySet(), subVolumeConductor_);
        problem_->bind(this->dipoleElement(), this->localDipolePosition(), this->dipole().moment());
        lop_ = std::make_shared<LOP>(*problem_, weighting_, config_.get<unsigned int>("intorderadd_eeg_patch"),
                                     config_.get<unsigned int>("intorderadd_eeg_boundary"));
        subFS_ = std::make_shared<SUBFS>(subVolumeConductor_);
        dataTree.set("sub_dofs", subFS_->getGFS().size());
        x_ = std::make_shared<DOF>(subFS_->getGFS(), 0.0);
        r_ = std::make_shared<DOF>(subFS_->getGFS(), 0.0);
        assembler_ = std::make_shared<AS>(*subFS_, *lop_, 1);
        timer.lap("sub_problem");
        timer.stop("bind_accumulated");
        // note: maybe invert normal in boundary condition of subtraction operator??
      }
      else if(continuityType == ContinuityType::continuous)
      {
        // create chi function
        HostLFS lfs(functionSpace_->getGFS());
        HostLFSCache indexMapper(lfs);

        // we first create the underlying vector containing the coefficients of chi when written in the FEM basis
        chiBasisCoefficientsPtr_ = std::make_shared<DOFVector>(functionSpace_->getGFS(), 0.0);

        // we now iterate over the inner region and set all coefficients corresponding to DOFs of elements contained in the inner region to 1.0
        for(const auto& element : patchAssembler_.patchElements()) {
          // bind local finite element space
          lfs.bind(element);
          indexMapper.update();

          for(size_t i = 0; i < indexMapper.size(); ++i) {
            auto index = indexMapper.containerIndex(i);
            (*chiBasisCoefficientsPtr_)[index] = 1.0;
          } // end inner for loop
        } //end outer for loop

        // we can now wrap chi into a grid function
        chiFunctionPtr_ = std::make_shared<DiscreteGridFunction>(functionSpace_->getGFS(), *chiBasisCoefficientsPtr_);
      } // end if
    } // end bind

    template<class Index>
    void print_entry(const std::pair<const Index, double>& index_value_pair, size_t& counter, double& current_max) const {
      double current_val = index_value_pair.second;
      if(std::abs(current_max) < std::abs(current_val)) {
        current_max = current_val;
      }
      std::cout << "val " << counter << " : \t" << index_value_pair.second << "\n";
      ++counter;
    }
    void print_entry(const double& value, size_t& counter, double& current_max) const {
      if(std::abs(current_max) < std::abs(value)) {
        current_max = value;
      }
      if(value != 0.0) {
        std::cout << "val " << counter << " : \t" << value << "\n";
      }
      ++counter;
    }

    virtual void assembleRightHandSide(VectorType& vector) const override
    {
      if constexpr(continuityType == ContinuityType::discontinuous)
      {
        assembleLocalDefaultSubtraction(vector);
        using LOP = LocalizedSubtractionDGLocalOperator<HostProblem, EdgeNormProvider, PenaltyFluxWeighting>;
        LOP lop2(*hostProblem_, edgeNormProvider_, weighting_, penalty_, intorderadd_eeg_boundary_);
        patchAssembler_.assemblePatchBoundary(vector, lop2);
      }
      else if(continuityType == ContinuityType::continuous)
      {
        using LOP = typename std::conditional<isP1FEM<FS>::value && dim == 3,
                                              LocalizedSubtractionCGP1LocalOperator<VC, DiscreteGridFunction, HostProblem>,
                                              LocalizedSubtractionCGLocalOperator<VC, DiscreteGridFunction, HostProblem>>::type;
        LOP cg_local_operator(volumeConductor_, chiFunctionPtr_, *hostProblem_, intorderadd_eeg_patch_, intorderadd_eeg_boundary_, intorderadd_eeg_transition_);
        patchAssembler_.assemblePatchVolume(vector, cg_local_operator);
        patchAssembler_.assemblePatchBoundary(vector, cg_local_operator);
        patchAssembler_.assembleTransitionVolume(vector, cg_local_operator);

        /*
        std::cout << " Entries of vector : \n";
        size_t counter = 0;
        double current_max = -1.0;
        for(const auto& entry : Dune::PDELab::Backend::native(vector)) {
          print_entry(entry, counter, current_max);
        }
        std::cout << " Largest value : \t" << current_max << "\n";
        */
      }
    }

    // NOTE : If you get weird results when enabling postprocessing for the localized subtraction source model,
    //        this might be a consequence of enabling the "subtract_mean" flag. If the dipole happens to lie exactly on
    //        or very close to a vertex, this value becomes NaN, or something very large. While the first directly
    //        shows you that something went wrong, the second can be more insideous, as adding a very large value simply deletes
    //        information in the smaller summand. Hence if your results do not match your expections, you might try to
    //        disable this flag and see if that helps. Note that this only applies to the direct forward solution and not
    //        to the transfer matrix approach. Hence if you use a direct approach I recommend doing rereferencing yourself
    //        after reducing to the electrode potentials, or whatever potentials you are interested in.
    virtual void postProcessSolution(VectorType& vector) const override
    {
      if constexpr(continuityType == ContinuityType::discontinuous)
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
      else if(continuityType == ContinuityType::continuous)
      {
        HostLFS lfs(functionSpace_->getGFS());
        HostLFSCache indexMapper(lfs);

        // wrap u infinity in GridViewFunction, enabling local interpolation
        auto u_inf = [this] (const CoordinateType& globalCoord) -> CoordinateField {return hostProblem_->get_u_infty(globalCoord);};
        auto gridFunctionUInfinity = Dune::Functions::makeGridViewFunction(u_inf, volumeConductor_->gridView());
        auto localUInfinity = localFunction(gridFunctionUInfinity);

        std::unordered_set<typename HostLFSCache::ContainerIndex> visited_dofs;
        for(const auto& element : patchAssembler_.patchElements()) {
          localUInfinity.bind(element);
          lfs.bind(element);
          indexMapper.update();
          std::vector<CoordinateField> u_infinity_interpolation(indexMapper.size());
          lfs.finiteElement().localInterpolation().interpolate(localUInfinity, u_infinity_interpolation);

          for(size_t i = 0; i < indexMapper.size(); ++i) {
            auto container_index = indexMapper.containerIndex(i);
            if(visited_dofs.count(container_index) == 0) {
              visited_dofs.insert(container_index);
              /*
              std::cout << "\n U_corr val : \t" << vector[container_index] << "\n";
              std::cout << " U_inf val : \t" << u_infinity_interpolation[i] << "\n";
              std::cout << " Sum : \t" <<  vector[container_index] + u_infinity_interpolation[i] << "\n\n";
              */
              vector[container_index] += u_infinity_interpolation[i];
            }
          }
        } // end loop over patch elements
      }
    } // end postProcessSolution

    virtual void
    postProcessSolution(const std::vector<CoordinateType>& electrodes,
                        std::vector<typename VectorType::field_type>& vector) const override
    {
      // note: need to check if electrode is within the patch
      // currently assume that the patch does not touch the boundary
    }

    virtual void postProcessMEG(const std::vector<CoordinateType>& coils,
                                const std::vector<std::vector<CoordinateType>>& projections,
                                std::vector<typename V::field_type>& fluxes) const
    {
      if constexpr(continuityType == ContinuityType::continuous) {
        fluxFromPatch(coils, projections, fluxes);
        fluxFromPatchBoundary(coils, projections, fluxes);
        fluxFromTransition(coils, projections, fluxes);
      }
      else {
        std::cout << " Noop postprocess\n";
      }
    }

  private:
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<const FS> functionSpace_;
    std::shared_ptr<SubVolumeConductor> subVolumeConductor_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<HostProblem> hostProblem_;
    std::shared_ptr<LOP> lop_;
    std::shared_ptr<SUBFS> subFS_;
    std::shared_ptr<AS> assembler_;
    std::shared_ptr<DOF> x_;
    std::shared_ptr<DOF> r_;
    Dune::ParameterTree config_;
    ElementPatchAssembler<VC, FS> patchAssembler_;
    EdgeNormProvider edgeNormProvider_;
    PenaltyFluxWeighting weighting_;
    unsigned int intorderadd_eeg_patch_;
    unsigned int intorderadd_eeg_boundary_;
    unsigned int intorderadd_eeg_transition_;
    unsigned int intorder_meg_patch_;
    unsigned int intorder_meg_boundary_;
    unsigned int intorder_meg_transition_;
    double penalty_;
    std::shared_ptr<DOFVector> chiBasisCoefficientsPtr_;
    std::shared_ptr<DiscreteGridFunction> chiFunctionPtr_;
    
    bool useAnalyticRHS_;

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

    // For the patch integrals we cannot assume the ratio of dipole distance and the edge length to be comfortably bounded below.
    // For low resolution meshes and highly eccentric sources this ratio might get very small. To compute the corresponding integrals
    // accurately we need appropriately high quadrature rules. The following function uses the largest edge length of an element to
    // choose such an order. The values should suffice under the following constraints:
    //    - the edge lengths of the elements do not exceed 6mm
    //    - the sources are at least 1mm from a conductivity jump
    // Note however that the integration orders were chosen using very conservative estimates, so you are probably fine even if your sources
    // are a little closer than 1mm to the conductivity jump, or the edge length exceeds 6mm by a little
    size_t intorderFromGeometry(const typename ElementType::Geometry& element_geometry) const
    {
      CoordinateField max_edge_length = -1.0;
      for(size_t i = 0; i < element_geometry.corners(); ++i) {
        for(size_t j = i + 1; j < element_geometry.corners(); ++j) {
          CoordinateField edge_length = (element_geometry.corner(i) - element_geometry.corner(j)).two_norm();
          if(edge_length > max_edge_length) max_edge_length = edge_length;
        }
      }

      if(max_edge_length <= 2.0) {
        return 8;
      }
      else if (max_edge_length <= 2.5) {
        return 9;
      }
      else if (max_edge_length <= 3) {
        return 11;
      }
      else if (max_edge_length <= 4) {
        return 13;
      }
      else {
        return 20;
      }
    }

    //////////////////////////////////////////////////
    // compute integral (sigma_corr grad(u_infinity)) x (x - y) / |x - y|^3 dy over patch region
    //////////////////////////////////////////////////
    void fluxFromPatch(const std::vector<CoordinateType>& coils,
                       const std::vector<std::vector<CoordinateType>>& projections,
                       std::vector<typename V::field_type>& fluxes) const
    {
      using GradientType = typename InfinityPotentialGradient<typename VC::GridView, CoordinateField>::RangeType;
      Tensor sigma_infinity = hostProblem_->get_sigma_infty();

      // loop over all patch elements and assemble local integrals
      for(const auto& element : patchAssembler_.patchElements()) {
        Tensor sigma = volumeConductor_->tensor(element);
        // elements with sigma_corr == 0 can be skipped
        if(sigma == sigma_infinity) {
          continue;
        }
        Tensor sigma_corr = sigma;
        sigma_corr -= sigma_infinity;

        const auto& elem_geo = element.geometry();
        auto local_coords_dummy = referenceElement(elem_geo).position(0, 0);
        auto jacobian_determinant = elem_geo.integrationElement(local_coords_dummy);

        // choose quadrature rule
        size_t intorder = (intorder_meg_patch_ == 0) ? intorderFromGeometry(elem_geo) : intorder_meg_patch_;
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim>& quad_rule = Dune::QuadratureRules<CoordinateField, dim>::rule(geo_type, intorder);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = jacobian_determinant * quad_point.weight();
          auto local_position = quad_point.position();
          auto global_position = elem_geo.global(local_position);

          // compute sigma_corr grad(u_infinity) * weight
          GradientType grad_u_infinity = hostProblem_->get_grad_u_infty(global_position);
          GradientType sigma_corr_grad_u_infinity;
          sigma_corr.mv(grad_u_infinity, sigma_corr_grad_u_infinity);
          sigma_corr_grad_u_infinity *= integration_factor;

          // loop over all coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            // compute (x - y) / |x - y|^3
            CoordinateType rhs = coils[i];
            rhs -= global_position;
            auto diff_norm = rhs.two_norm();
            auto norm_cubed = diff_norm * diff_norm * diff_norm;
            rhs /= norm_cubed;

            // compute cross product
            CoordinateType cross_product;
            Dune::PDELab::CrossProduct<dim, dim>(cross_product, sigma_corr_grad_u_infinity, rhs);
            for(size_t j = 0; j < projections[i].size(); ++j) {
              fluxes[counter] += cross_product * projections[i][j];
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quadrature points
      } // end loop over patch elements
    } // end fluxFromPatch

    //////////////////////////////////////////////////
    // compute integral (sigma grad(chi * u_infinity)) x (x - y) / |x - y|^3 dy over transition region
    //////////////////////////////////////////////////
    void fluxFromTransition(const std::vector<CoordinateType>& coils,
                            const std::vector<std::vector<CoordinateType>>& projections,
                            std::vector<typename V::field_type>& fluxes) const
    {
      using GradientType = typename InfinityPotentialGradient<typename VC::GridView, CoordinateField>::RangeType;

      LocalFunction chi_local = localFunction(*chiFunctionPtr_);
      LocalDerivativeFunction grad_chi_local = localFunction(derivative(*chiFunctionPtr_));

      for(const auto& element : patchAssembler_.transitionElements()) {
        Tensor sigma = volumeConductor_->tensor(element);

        chi_local.bind(element);
        grad_chi_local.bind(element);

        const auto& elem_geo = element.geometry();
        auto local_coords_dummy = referenceElement(elem_geo).position(0, 0);
        auto jacobian_determinant = elem_geo.integrationElement(local_coords_dummy);

        // choose quadrature rule
        Dune::GeometryType geo_type = elem_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim>& quad_rule = Dune::QuadratureRules<CoordinateField, dim>::rule(geo_type, intorder_meg_transition_);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = jacobian_determinant * quad_point.weight();
          auto local_position = quad_point.position();
          auto global_position = elem_geo.global(local_position);

          // compute u_infinity * grad_chi
          GradientType u_infinity_grad_chi = grad_chi_local(local_position)[0];
          u_infinity_grad_chi *= hostProblem_->get_u_infty(global_position);

          // compute chi * grad_u_infinity
          GradientType chi_grad_u_infinity = hostProblem_->get_grad_u_infty(global_position);
          chi_grad_u_infinity *= chi_local(local_position);

          // compute LHS of cross product
          u_infinity_grad_chi += chi_grad_u_infinity;
          GradientType lhs;
          sigma.mv(u_infinity_grad_chi, lhs);
          lhs *= integration_factor;

          // loop over all coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            for(size_t j = 0; j < projections[i].size(); ++j) {
              // compute RHS
              CoordinateType rhs = coils[i] - global_position;
              auto diff_norm = rhs.two_norm();
              auto norm_cubed = diff_norm * diff_norm * diff_norm;
              rhs /= norm_cubed;

              // compute cross product
              GradientType crossProduct;
              Dune::PDELab::CrossProduct<dim, dim>(crossProduct, lhs, rhs);

              fluxes[counter] += (crossProduct * projections[i][j]);
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quadrature points
      } // end loop over transition elements
    } // end fluxFromTransition


    //////////////////////////////////////////////////
    // compute integral sigma_infinity u_infinity (eta x (x - y)/ |x - y|^3) ds
    //////////////////////////////////////////////////
    void fluxFromPatchBoundary(const std::vector<CoordinateType>& coils,
                               const std::vector<std::vector<CoordinateType>>& projections,
                               std::vector<typename V::field_type>& fluxes) const
    {
      // iterate over all boundary patch boundary intersections
      for(const auto& intersection : patchAssembler_.intersections()) {
        // get intersection geometry
        const auto& intersection_geo = intersection.geometry();
        auto local_coords_dummy = referenceElement(intersection_geo).position(0, 0);
        auto gramian = intersection_geo.integrationElement(local_coords_dummy);
        CoordinateType unitOuterNormal = intersection.centerUnitOuterNormal();

        Tensor sigma_infinity = hostProblem_->get_sigma_infty();

        // choose quadrature rule
        Dune::GeometryType geo_type = intersection_geo.type();
        const Dune::QuadratureRule<CoordinateField, dim - 1>& quad_rule = Dune::QuadratureRules<CoordinateField, dim - 1>::rule(geo_type, intorder_meg_boundary_);

        // perform the integration
        for(const auto& quad_point : quad_rule) {
          auto integration_factor = gramian * quad_point.weight();
          auto local_position = quad_point.position();
          auto global_position = intersection_geo.global(local_position);

          // compute sigma_infinity_u_infinity
          // NOTE : assumes isotropic sigma_infinity
          auto sigma_infinity_u_infinity = sigma_infinity[0][0] * hostProblem_->get_u_infty(global_position);
          sigma_infinity_u_infinity *= integration_factor;

          // loop over coils and projections
          size_t counter = 0;
          for(size_t i = 0; i < coils.size(); ++i) {
            for(size_t j = 0; j < projections[i].size(); ++j) {
              // compute cross product
              CoordinateType rhs = coils[i] - global_position;
              auto diff_norm = rhs.two_norm();
              auto norm_cubed = diff_norm * diff_norm * diff_norm;
              rhs /= norm_cubed;
              CoordinateType crossProduct;
              Dune::PDELab::CrossProduct<dim, dim>(crossProduct, unitOuterNormal, rhs);

              fluxes[counter] += sigma_infinity_u_infinity * (crossProduct * projections[i][j]);
              ++counter;
            } // end loop over projections
          } // end loop over coils
        } // end loop over quadrature points
      } // end loop over intersections
    } // end fluxFromPatchBoundary

  };
}

#endif // DUNEURO_LOCALIZED_SUBTRACTION_SOURCE_MODEL_HH
