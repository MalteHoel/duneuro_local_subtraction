// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_IO_VOLUME_CONDUCTOR_VTK_WRITER_HH
#define DUNEURO_IO_VOLUME_CONDUCTOR_VTK_WRITER_HH

/*
 * This class implements a writer for a volume conductor. It should do the following:
 *    -) Print the volume conductor, which means the underlying mesh and the associated conductivities
 *    -) It should be able to visualize FEM trial functions on the volume conductor
 */

#include <string>

#include <duneuro/common/function.hh>
#include <duneuro/io/vtk_writer.hh>
#include <duneuro/io/fitted_tensor_vtk_functor.hh>

#if HAVE_DUNE_UDG
#include <duneuro/io/refined_vtk_writer.hh>
#endif

#include <duneuro/io/vtk_functors.hh>

namespace duneuro {

  class VolumeConductorVTKWriterInterface {
  public:
    virtual void addVertexData(const Function& function, const std::string& name) = 0;
    virtual void addVertexDataGradient(const Function& function, const std::string& name) = 0;
    virtual void addCellData(const Function& function, const std::string& name) = 0;
    virtual void addCellDataGradient(const Function& function, const std::string& name) = 0;
    virtual void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) = 0;
    
    virtual ~VolumeConductorVTKWriterInterface() {}
  };

  template<class Solver>
  class VolumeConductorVTKWriter
    : public VolumeConductorVTKWriterInterface   
  {
  public:
    using VolumeConductor = typename Solver::Traits::VolumeConductor;
    enum {dim = VolumeConductor::dim};
    using DOFVector = typename Solver::Traits::DomainDOFVector;
    using GridFunctionSpace = typename Solver::Traits::FunctionSpace::GFS;
    using DiscreteGridFunction = typename Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, DOFVector>;
    using VTKFunctionAdapter = typename Dune::PDELab::VTKGridFunctionAdapter<DiscreteGridFunction>;
    using DGFGradient = typename Dune::PDELab::DiscreteGridFunctionGradient<GridFunctionSpace, DOFVector>;
    using VTKGradientAdapter = typename Dune::PDELab::VTKGridFunctionAdapter<DGFGradient>;
    using Writer = VTKWriter<typename Solver::Traits::GridView>;
  
    VolumeConductorVTKWriter(const Solver& solver, bool visualizeAnisotropy)
      : solver_(solver)
      , gfs_(solver_.functionSpace().getGFS())
      , writer_(solver.volumeConductor()->gridView())
    {
      writer_.addCellData(std::make_shared<FittedLabelFunctor<VolumeConductor>>(solver_.volumeConductor()));
      writer_.addCellData(std::make_shared<FittedTensorNormFunctor<VolumeConductor>>(solver_.volumeConductor()));
      if(visualizeAnisotropy) {
        writer_.addCellData(std::make_shared<FittedTensorFunctor<VolumeConductor>>(solver_.volumeConductor()));
        writer_.addCellData(std::make_shared<FittedTensorFractionalAnisotropyFunctor<VolumeConductor>>(solver_.volumeConductor()));
#if HAVE_EIGEN
        for (unsigned int i = 0; i < dim; ++i) {
          writer_.addCellData(std::make_shared<FittedTensorEigenvectorFunctor<VolumeConductor>>(solver_.volumeConductor(), i));
        }
#endif
      }
    }
    
    virtual void addVertexData(const Function& function, const std::string& name) override
    {
      std::shared_ptr<DiscreteGridFunction> gridFunctionPtr = std::make_shared<DiscreteGridFunction>(gfs_, function.cast<DOFVector>());
      std::shared_ptr<VTKFunctionAdapter> vtkFunctionPtr = std::make_shared<VTKFunctionAdapter>(gridFunctionPtr, name);
      writer_.addVertexData(vtkFunctionPtr);
    }
    
    virtual void addVertexDataGradient(const Function& function, const std::string& name) override
    {
      std::shared_ptr<DGFGradient> gridFunctionPtr = std::make_shared<DGFGradient>(gfs_, function.cast<DOFVector>());
      std::shared_ptr<VTKGradientAdapter> vtkFunctionPtr = std::make_shared<VTKGradientAdapter>(gridFunctionPtr, name);
      writer_.addVertexData(vtkFunctionPtr);
    }
    
    virtual void addCellData(const Function& function, const std::string& name) override
    {
      std::shared_ptr<DiscreteGridFunction> gridFunctionPtr = std::make_shared<DiscreteGridFunction>(gfs_, function.cast<DOFVector>());
      std::shared_ptr<VTKFunctionAdapter> vtkFunctionPtr = std::make_shared<VTKFunctionAdapter>(gridFunctionPtr, name);
      writer_.addCellData(vtkFunctionPtr);
    }
    
    virtual void addCellDataGradient(const Function& function, const std::string& name) override
    {
      std::shared_ptr<DGFGradient> gridFunctionPtr = std::make_shared<DGFGradient>(gfs_, function.cast<DOFVector>());
      std::shared_ptr<VTKGradientAdapter> vtkFunctionPtr = std::make_shared<VTKGradientAdapter>(gridFunctionPtr, name);
      writer_.addCellData(vtkFunctionPtr);
    }
    
    virtual void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree()) override
    {
      writer_.write(config, dataTree);
    }
    
  private:
    const Solver& solver_;
    const GridFunctionSpace& gfs_;
    Writer writer_;
  };


#if HAVE_DUNE_UDG
  template<class Solver>
  class UnfittedVCVTKWriter
    : public VolumeConductorVTKWriterInterface
  {
  public:
    using GridView = typename Solver::Traits::GridView;
    using DOFVector = typename Solver::Traits::DomainDOFVector;
    using GridFunctionSpace = typename Solver::Traits::FunctionSpace::GFS;
    using SubTriangulation = typename Solver::Traits::SubTriangulation;
    enum {compartments = Solver::Traits::compartments};
    using Writer = RefinedVTKWriter<GridFunctionSpace, SubTriangulation, compartments>;
  
    using TensorVTKFunction = TensorUnfittedVTKGridFunction<GridView>;
    using DomainIndexVTKFunction = typename Dune::UDG::DomainIndexUnfittedVTKGridFunction<GridView>;
    using HostIndexVTKFunction = typename Dune::UDG::HostCellIndexUnfittedVTKGridFunction<GridView>;
  
    UnfittedVCVTKWriter(const std::shared_ptr<const Solver> solverPtr,  std::shared_ptr<SubTriangulation> subTriangulationPtr, const GridView& fundamentalGridView, std::vector<double> conductivities, const std::string& modeString, bool scaleToBBox = true)
      : solverPtr_(solverPtr)
      , fundamentalGridView_(fundamentalGridView)
      , writer_(subTriangulationPtr, solverPtr_->functionSpace().getGFS(), scaleToBBox)
    {
      writer_.addVertexData(std::make_shared<TensorVTKFunction>(fundamentalGridView_, conductivities));
      writer_.addVertexData(std::make_shared<DomainIndexVTKFunction>(fundamentalGridView_));
      
      if((modeString == "faces") || (modeString == "boundary")) {
        writer_.addVertexData(std::make_shared<DomainIndexVTKFunction>(fundamentalGridView_, false));
      }
      else {
        writer_.addVertexData(std::make_shared<HostIndexVTKFunction>(fundamentalGridView_));
      }
    }
    
    virtual void addVertexData(const Function& function, const std::string& name)
    {
      writer_.addVertexData(*solverPtr_, function.cast<DOFVector>(), name);
    }
    
    virtual void addVertexDataGradient(const Function& function, const std::string& name)
    {
      writer_.addVertexDataGradient(*solverPtr_, function.cast<DOFVector>(), name);
    }
    
    virtual void addCellData(const Function& function, const std::string& name)
    {
      writer_.addCellData(*solverPtr_, function.cast<DOFVector>(), name);
    }
    
    virtual void addCellDataGradient(const Function& function, const std::string& name)
    {
      DUNE_THROW(Dune::Exception, "adding cell data gradient is currently not supported for unfitted volume conductors, use vertex data gradient instead");
    }
    
    virtual void write(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      writer_.write(config, dataTree);
    }
   
  private:
    const std::shared_ptr<const Solver> solverPtr_;
    const GridView& fundamentalGridView_;
    Writer writer_;
   };
#endif

} // namespace duneuro
#endif // DUNEURO_IO_VOLUME_CONDUCTOR_VTK_WRITER_HH
