#ifndef DUNEURO_CUTFEM_TDCS_DRIVER_HH
#define DUNEURO_CUTFEM_TDCS_DRIVER_HH

#include <dune/common/parametertree.hh>

#include <dune/udg/simpletpmctriangulation.hh>

#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/structured_grid_utilities.hh>
#include <duneuro/tes/tdcs_driver_interface.hh>
#include <duneuro/tes/tdcs_patch_udg_parameter.hh>
#include <duneuro/tes/tdcs_rhs_factory.hh>
#include <duneuro/tes/tdcs_solver.hh>

#if HAVE_TBB
#include <tbb/tbb.h>
#endif
#include <dune/udg/pdelab/gridfunction.hh>
#include <dune/common/std/memory.hh>
#include <dune/common/version.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelement.hh>

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

#include <duneuro/common/cutfem_solver.hh>
#include <duneuro/common/cutfem_solver_backend.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/stl.hh>
#include <duneuro/eeg/projected_electrodes.hh>
#include <duneuro/io/refined_vtk_writer.hh>
#include <duneuro/io/vtk_functors.hh>
#include <duneuro/meeg/unfitted_meeg_driver_data.hh>
#include <duneuro/udg/subtriangulation_statistics.hh>


namespace duneuro
{
  template <int dim, int degree, int compartments>
  struct CutFEMTDCSDriverTraits {
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GridView = typename Grid::LevelGridView;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using Problem = TDCSPatchUDGParameter<GridView>; 
    using Solver = CutFEMSolver<SubTriangulation, compartments, degree, Problem>;
    using SolverBackend = CutFEMSolverBackend<Solver>;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
  };

  template <int dim, int degree, int compartments>
  class CutFEMTDCSDriver : public TDCSDriverInterface<dim>
  {
  public:
    using Traits = CutFEMTDCSDriverTraits<dim, degree, compartments>;

    explicit CutFEMTDCSDriver(const PatchSet<double, dim>& patchSet, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : CutFEMTDCSDriver(UnfittedMEEGDriverData<dim>{}, patchSet, config, dataTree)
    {
    }

    explicit CutFEMTDCSDriver(const UnfittedMEEGDriverData<dim>& data,
                           const PatchSet<double, dim>& patchSet, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : config_(config)
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false),
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 6)
              Dune::UDG::simpleTPMCIntersectionsFromString(
                  config.get<std::string>("udg.intersections", "all")),
#endif
              config.get<double>("udg.value_tolerance", 1e-8)))
        , problem_(std::make_shared<typename Traits::Problem>(
              config_.get<std::vector<double>>("solver.conductivities"), patchSet))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(fundamentalGridView_))
        , solver_(std::make_shared<typename Traits::Solver>(subTriangulation_, elementSearch_,
                                                            problem_, config.sub("solver")))
        , solverBackend_(std::make_shared<typename Traits::SolverBackend>(
              solver_, config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree()))
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
    {
    }

    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }
    virtual std::unique_ptr<DenseMatrix<double>> CenterEvaluation(const Function& solution)
    {}

    virtual void solveTDCSForward(Function& solution, const Dune::ParameterTree& config,
                                  DataTree dataTree = DataTree()) override
    {
      solver_->solve(solverBackend_->get(), solution.cast<typename Traits::DomainDOFVector>(),
                     config, dataTree);
    }

    virtual void write(const Function& function, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS());
        vtkWriter.addVertexData(*solver_, function.cast<typename Traits::DomainDOFVector>(),
                                "potential");
        vtkWriter.addVertexDataGradient(*solver_, function.cast<typename Traits::DomainDOFVector>(),
                                        "gradient_potential");
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        auto modeString = config.get<std::string>("mode", "volume");
        if ((modeString == "faces") || (modeString == "boundary")) {
          vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                      typename Traits::GridView>>(fundamentalGridView_, false));
        }
        vtkWriter.write(config, dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual void write(const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS());
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        auto modeString = config.get<std::string>("mode", "volume");
        if ((modeString == "faces") || (modeString == "boundary")) {
          vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                      typename Traits::GridView>>(fundamentalGridView_, false));
        }
        vtkWriter.write(config, dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }
    virtual void
    setElectrodes(const std::vector<typename TDCSDriverInterface<dim>::CoordinateType>& electrodes,
                  const Dune::ParameterTree& config) override
    {}
    virtual std::unique_ptr<DenseMatrix<double>>
    computeEvaluationMatrix(const Dune::ParameterTree& config,
                             DataTree dataTree = DataTree()) override{}
  
    virtual std::vector<std::vector<double>> applyEvaluationMatrix(const DenseMatrix<double>& EvaluationMatrix,
                                           const std::vector<typename TDCSDriverInterface<dim>::CoordinateType>& positions,
                                           Dune::ParameterTree cfg, DataTree dataTree = DataTree() ) const override {}


    virtual std::vector<std::vector<double>> applyEvaluationMatrix(const DenseMatrix<double>& EvaluationMatrix,
                                           Dune::ParameterTree cfg, DataTree dataTree = DataTree() ) const override {}
  private:
    Dune::ParameterTree config_;
    std::unique_ptr<typename Traits::Grid> grid_;
    typename Traits::GridView fundamentalGridView_;
    typename Traits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView> domain_;
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Problem> problem_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<typename Traits::SolverBackend> solverBackend_;
    std::vector<double> conductivities_;
  };


    //Cutfem-Driver for Point-Electrodes


  template <int dim, int degree, int compartments>
  struct CutFEMTDCSPointDriverTraits {
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GridView = typename Grid::LevelGridView;
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GridView, GridView>;
    using ElementSearch = KDTreeElementSearch<GridView>;
    using Solver = CutFEMSolver<SubTriangulation, compartments, degree>;
    using SolverBackend = CutFEMSolverBackend<Solver>;
    using ProjectedPosition = duneuro::ProjectedElectrode<GridView>;
    using DomainDOFVector = typename Solver::Traits::DomainDOFVector;
    using RangeDOFVector = typename Solver::Traits::RangeDOFVector;
    using TdcsRHSFactory = UnfittedTdcsRHSFactory;
  };

  template <int dim, int degree, int compartments>
  class CutFEMTDCSPointDriver : public TDCSDriverInterface<dim>
  {
  public:
    using Traits = CutFEMTDCSPointDriverTraits<dim, degree, compartments>;

    explicit CutFEMTDCSPointDriver( const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : CutFEMTDCSPointDriver(UnfittedMEEGDriverData<dim>{}, config, dataTree)
    {
    }

    explicit CutFEMTDCSPointDriver(const UnfittedMEEGDriverData<dim>& data,
                            const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : data_(data)
        , config_(config)
        , grid_(make_structured_grid<dim>(config.sub("volume_conductor.grid")))
        , fundamentalGridView_(grid_->levelGridView(0))
        , levelSetGridView_(grid_->levelGridView(grid_->maxLevel()))
        , domain_(levelSetGridView_, data_.levelSetData, config.sub("domain"))
        , subTriangulation_(std::make_shared<typename Traits::SubTriangulation>(
              fundamentalGridView_, levelSetGridView_, domain_.getDomainConfiguration(),
              config.get<bool>("udg.force_refinement", false),
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 6)
              Dune::UDG::simpleTPMCIntersectionsFromString(
                  config.get<std::string>("udg.intersections", "all")),
#endif
              config.get<double>("udg.value_tolerance", 1e-8)))
        , elementSearch_(std::make_shared<typename Traits::ElementSearch>(fundamentalGridView_))
        , solver_(std::make_shared<typename Traits::Solver>(subTriangulation_, elementSearch_,
                                                             config.sub("solver")))
        , conductivities_(config.get<std::vector<double>>("solver.conductivities"))
        , solverBackend_(solver_,
                         config.hasSub("solver") ? config.sub("solver") : Dune::ParameterTree())     
        , rightHandSideVector_(make_range_dof_vector(*solver_, 0.0))
        , tdcsSolver_(solver_, *subTriangulation_, config_)
    {
    }

    virtual std::unique_ptr<Function> makeDomainFunction() const override
    {
      return Dune::Std::make_unique<Function>(make_domain_dof_vector(*solver_, 0.0));
    }
  

    virtual void
    setElectrodes(const std::vector<typename TDCSDriverInterface<dim>::CoordinateType>& electrodes,
                  const Dune::ParameterTree& config) override
    {
      projectedElectrodes_ = Dune::Std::make_unique<ProjectedElectrodes<typename Traits::GridView>>(
          electrodes, solver_->functionSpace().getGFS(), *subTriangulation_);
      projectedGlobalElectrodes_.clear();
      for (unsigned int i = 0; i < projectedElectrodes_->size(); ++i) {
        projectedGlobalElectrodes_.push_back(projectedElectrodes_->projection(i));
      }
    }

    virtual void solveTDCSForward(Function& solution, const Dune::ParameterTree& config,
                                  DataTree dataTree = DataTree()) override
    {
      
      
      // assemble right hand side
      Dune::Timer timer;
      auto Ref = projectedElectrodes_->getProjection(0);
      auto Elec = projectedElectrodes_->getProjection(1);
      auto rhsAssembler =
      UnfittedTdcsRHSFactory::template create<typename Traits::RangeDOFVector>(*solver_, *subTriangulation_, config);  
      rhsAssembler->bind(projectedElectrodes_->element(0), Ref.localPosition, projectedElectrodes_->element(1),
                         Elec.localPosition);
      rhsAssembler->assembleRightHandSide(*rightHandSideVector_);          
      timer.stop();
      dataTree.set("time_rhs_assembly", timer.lastElapsed());
      timer.start();

      // solve system
      solver_->solve(solverBackend_.get(),*rightHandSideVector_, solution.cast<typename Traits::DomainDOFVector>(), config.sub("solver"),
                     dataTree.sub("linear_system_solver"));
      timer.stop();
      dataTree.set("time_solve", timer.lastElapsed());
      dataTree.set("time", timer.elapsed());

    }


    virtual std::unique_ptr<DenseMatrix<double>>
    computeEvaluationMatrix(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) override
    {
      return tdcsSolver_.tdcsEvaluationMatrix(solverBackend_, *projectedElectrodes_, config,
                                            dataTree);
    }
  
virtual std::vector<std::vector<double>> applyEvaluationMatrix(const DenseMatrix<double>& EvaluationMatrix,
                                           const std::vector<typename TDCSDriverInterface<dim>::CoordinateType>& positions,
                                           Dune::ParameterTree cfg, DataTree dataTree = DataTree() ) const override
{
  KDTreeElementSearch<typename Traits::GridView> search(solver_->functionSpace().getGFS().gridView());
  std::vector<typename TDCSDriverInterface<dim>::CoordinateType> localPositions(positions.size());
  std::vector<typename Traits::GridView::template Codim<0>::Entity> elements(positions.size());
  std::size_t index = 0;
  for (const auto& coord : positions)
  {
    const auto& element = search.findEntity(coord);
    elements[index] = element;
    const auto local = element.geometry().local(coord);
    localPositions[index] = local;
    index++;
  }

  return tdcsSolver_.applyEvaluationMatrix(EvaluationMatrix, elements, localPositions, cfg);
}


virtual std::vector<std::vector<double>> applyEvaluationMatrix(const DenseMatrix<double>& EvaluationMatrix,
                                           Dune::ParameterTree cfg, DataTree dataTree = DataTree() ) const override
{
  unsigned int numberHostCells = 0;
    for (const auto& element : Dune::elements(fundamentalGridView_))    // Determine number of Host Cells
    {
        if (subTriangulation_->isHostCell(element)) {
          numberHostCells+=1;
        }
    }
  std::vector<typename TDCSDriverInterface<dim>::CoordinateType> localPositions(numberHostCells);
  std::vector<typename Traits::GridView::template Codim<0>::Entity> elements(numberHostCells);
  std::size_t index = 0;

  for (const auto& element : Dune::elements(fundamentalGridView_))    
  {
    if (!subTriangulation_->isHostCell(element))                      // skip elements that are outside the brain
    {
      continue;
    }
    elements[index] = element;
    const auto local = element.geometry().local(element.geometry().center());
    localPositions[index] = local;
    index++;
  }

 return tdcsSolver_.applyEvaluationMatrix(EvaluationMatrix, elements, localPositions, cfg);
}
 // to be deleted once a more suitable evaluation interface is implemented
  virtual std::unique_ptr<DenseMatrix<double>> CenterEvaluation(const Function& solution)
  {
    auto testMat = Dune::Std::make_unique<DenseMatrix<double>>(
        1,
        solver_->functionSpace().getGFS().ordering().size()); 
      set_matrix_row(*testMat, 0, Dune::PDELab::Backend::native(solution.cast<typename Traits::DomainDOFVector>()));


    unsigned int NumberHostCells = 0;

    for (const auto& element : Dune::elements(subTriangulation_->gridView()))    // Determine number of Host Cells
    {
        if (subTriangulation_->isHostCell(element)) {
          NumberHostCells+=1;
        }
    }

  //  KDTreeElementSearch<typename Traits::GridView> search(solver_->functionSpace().getGFS().gridView());
  //  unsigned int numberHostCells = EvalPoints.size();
    auto elementCenter = Dune::Std::make_unique<DenseMatrix<double>>(   // Matrix that will later be filled with the Center Positions, Potentials and Gradients
    NumberHostCells,2*Traits::GridView::dimension + 1);

    std::size_t offset = 0;
    for (const auto& element : Dune::elements(subTriangulation_->gridView())) {

    //const auto& element = search.findEntity(coord);

      if (!subTriangulation_->isHostCell(element))                      // skip elements that are outside the brain
      {
      //    offset+=1;
          continue;
      }
      auto y=element.geometry().center();
      //auto y = coord;
      auto localPos = element.geometry().local(y);                      
      auto pot = evaluateatCoordinate(element, localPos, solution.cast<typename Traits::DomainDOFVector>()); 
      for(std::size_t i=0; i<Traits::GridView::dimension; ++i)
      {
        (*elementCenter)(offset,i) = y[i];    // xyz-coordinates of the center
      }
      for(std::size_t i = 0;i<4;i++)
      {
        (*elementCenter)(offset,i+3) = pot[i];    // potential and gradient values
      }
      offset+=1;
    }
    return elementCenter;
  }

   // Evaluation of the Potential at the local Coordinate of an Element x
   // 2 Versions: 

   // The first creates Subtriangulation of the element,
   // binds the LFS of the different domains to the element, evaluates them at the local coordinate, multiplies them with 
   // the corresponding coeff. from the DOF-solution Vector and returns the sum.

   //The second uses  GridfunctionSubspace and DiscreteGridfunction.
   template<typename Element, typename LocalCoordinate>
   const std::vector<double> evaluateatCoordinate(Element& element, LocalCoordinate& local, const typename Traits::DomainDOFVector& solution)
      {  
        
      using GFS = typename Traits::Solver::Traits::FunctionSpace::GFS;
      using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
      using UST = Dune::PDELab::UnfittedSubTriangulation<typename Traits::GridView>;
      using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
      using ChildLFS = typename ULFS::template Child<0>::Type;
      using FESwitch = Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using RangeType = typename BasisSwitch::Range;
      using RangeFieldType = typename BasisSwitch::RangeField;
      using Real = typename Traits::GridView::ctype;
     


      ULFS ulfs(solver_->functionSpace().getGFS());
      UCache ucache(ulfs);
      UST ust(subTriangulation_->gridView(), *subTriangulation_);
      auto global = element.geometry().global(local); 
      int comp = 0;
      for (int i = 0;i<dim;i++)
      {
        global[i]-=127;
      }

      double ecc = global.two_norm();
      global[0]-=2;
      double ecc2 = global.two_norm();
      if (ecc<86)
      {comp = 1;}
      if (ecc<80)
      {
          comp = 2;
        }
      if (ecc2<78)
      {comp = 3;}
      Dune::FieldVector<Real, dim> y;
      std::vector<RangeType> phi;                                 // storage for Ansatzfunction values
      std::vector<double> output(4);      
      ust.create(element);                                        // splitting the Element 
      for (const auto& ep : ust) 
      {
        if(ep.domainIndex()!= comp)
        {continue;}
          ChildLFS& childLfs(ulfs.child(ep.domainIndex() ) );     // chooses the correct Ansatzfunctionspace and binds it to the El.
          ulfs.bind(ep.subEntity(), true);
          ucache.update();
          if (childLfs.size() == 0)
          {
            continue;
          }
          /*
          FESwitch::basis(childLfs.finiteElement()).reset();
        
          phi.resize(childLfs.size());
          std::vector<Dune::FieldMatrix<Real, 1, dim>> gradphi(childLfs.size());
          FESwitch::basis(childLfs.finiteElement()).evaluateFunction(local, phi);         // Ansatzfct eval.

          FESwitch::basis(childLfs.finiteElement()).evaluateJacobian(local, gradphi);     // Gradients
          for (unsigned int i = 0; i < ucache.size(); ++i) {
            const double& coeff = solution[ucache.containerIndex(childLfs.localIndex(i))];
            output[0] += phi[i]*coeff;
            output[1] += gradphi[i][0][0]*coeff;      // gradphi is a vector of 1xDim Matrices,so
            output[2] += gradphi[i][0][1]*coeff;      // i is the number of the Ansatzfunction and 
            output[3] += gradphi[i][0][2]*coeff;      // the first 0 is the row index
          }  
          */
                 // get Jacobian of geometry
                 
        const auto JgeoIT = element.geometry().jacobianInverseTransposed(local);
        // get local Jacobians/gradients of the shape functions
        std::vector<Dune::FieldMatrix<Real, 1, dim>> J(childLfs.size());
        FESwitch::basis(childLfs.finiteElement()).reset();
        FESwitch::basis(childLfs.finiteElement()).evaluateJacobian(local, J);     // Gradients

  
        Dune::FieldVector<Real, dim> gradphi;
        y = 0;
        for(unsigned int i = 0; i < ucache.size(); ++i) {
          // compute global gradient of shape function i
          gradphi = 0;
          JgeoIT.umv(J[i][0], gradphi);
          const double& coeff = solution[ucache.containerIndex(childLfs.localIndex(i))];
          // sum up global gradients, weighting them with the appropriate coeff
          y.axpy(coeff, gradphi);


        }

          break;
          
      }
      for (int i = 0;i<dim; ++i)
      {
        output[i+1] = y[i];
      }
      return output;

    // Alternative potential computation
    
/*
      using SubGFS =
          Dune::PDELab::GridFunctionSubSpace<typename Traits::Solver::Traits::FunctionSpace::GFS,
                                             Dune::TypeTree::TreePath<0>>;      // Does this 0 represent Domain Index 0 aka the skin?
      SubGFS subGFS(solver_->functionSpace().getGFS());
      using DGF = Dune::PDELab::DiscreteGridFunction<SubGFS, typename Traits::DomainDOFVector>;
      DGF dgf(subGFS, solution);
      double pot;   
      dgf.evaluate(element, local, pot);
     
      std::vector<double> output(4);
      output[0] = pot;
      using DGFG = Dune::PDELab::DiscreteGridFunctionGradient<SubGFS,
                                                              typename Traits::DomainDOFVector>;
      DGFG dgfg(subGFS, solution);
      typename DGFG::Traits::RangeType y;
      dgfg.evaluate(element, local, y); 
      for(std::size_t i = 0; i<3; i++)  
      {
        output[i+1] = y[i];
      }
      return output;
      */
      }

    virtual void write(const Function& function, const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS());
        vtkWriter.addVertexData(*solver_, function.cast<typename Traits::DomainDOFVector>(),
                                "potential");
        vtkWriter.addVertexDataGradient(*solver_, function.cast<typename Traits::DomainDOFVector>(),
                                        "gradient_potential");
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        auto modeString = config.get<std::string>("mode", "volume");
        if ((modeString == "faces") || (modeString == "boundary")) {
          vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                      typename Traits::GridView>>(fundamentalGridView_, false));
        }
        vtkWriter.write(config, dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }

    virtual void write(const Dune::ParameterTree& config,
                       DataTree dataTree = DataTree()) const override
    {
      auto format = config.get<std::string>("format");
      if (format == "vtk") {
        RefinedVTKWriter<typename Traits::Solver::Traits::FunctionSpace::GFS,
                         typename Traits::SubTriangulation, compartments>
            vtkWriter(subTriangulation_, solver_->functionSpace().getGFS());
        vtkWriter.addVertexData(
            std::make_shared<TensorUnfittedVTKGridFunction<typename Traits::GridView>>(
                fundamentalGridView_, conductivities_));
        vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                    typename Traits::GridView>>(fundamentalGridView_));
        auto modeString = config.get<std::string>("mode", "volume");
        if ((modeString == "faces") || (modeString == "boundary")) {
          vtkWriter.addVertexData(std::make_shared<Dune::UDG::DomainIndexUnfittedVTKGridFunction<
                                      typename Traits::GridView>>(fundamentalGridView_, false));
        }
        vtkWriter.write(config, dataTree);
      } else {
        DUNE_THROW(Dune::Exception, "Unknown format \"" << format << "\"");
      }
    }
  private:
    UnfittedMEEGDriverData<dim> data_;
    Dune::ParameterTree config_;
    std::unique_ptr<typename Traits::Grid> grid_;
    typename Traits::GridView fundamentalGridView_;
    typename Traits::GridView levelSetGridView_;
    SimpleTPMCDomain<typename Traits::GridView, typename Traits::GridView> domain_;
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::ElementSearch> elementSearch_;
    std::shared_ptr<typename Traits::Solver> solver_;
    typename Traits::SolverBackend solverBackend_;
    std::vector<double> conductivities_;
    std::unique_ptr<typename Traits::RangeDOFVector> rightHandSideVector_;
    std::unique_ptr<ProjectedElectrodes<typename Traits::GridView>> projectedElectrodes_;
    std::vector<Dune::FieldVector<typename Traits::GridView::ctype, Traits::GridView::dimension>>
        projectedGlobalElectrodes_;
    TDCSSolver<typename Traits::Solver, typename Traits::TdcsRHSFactory> tdcsSolver_;
  };
}

#endif // DUNEURO_CUTFEM_TDCS_DRIVER_HH
