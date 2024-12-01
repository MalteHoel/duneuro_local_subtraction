// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_DOF_VECTOR_EVALUATOR_HH
#define DUNEURO_DOF_VECTOR_EVALUATOR_HH

#include <vector>
#include <memory>
#include <type_traits>
#include <string>
#include <limits>

#include <dune/common/shared_ptr.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include <duneuro/common/dense_matrix.hh>
#include <duneuro/common/matrix_utilities.hh>

#if HAVE_DUNE_UDG
#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>
#include <dune/udg/pdelab/gridfunction.hh>
#endif

namespace duneuro {

  /* 
   * This class is supposed to aid in the evaluation of DOF vectors, given in various forms.
   * A DOF vector defines a function on the underlying grid. Depending on the question at hand, 
   * one might now want to evaluation this function, or its derivative, at some predefined positions.
   * This class implements this evaluation. We aim to support the evaluation of the function itself,
   * its gradient, and - sigma * gradient of the function (which typically corresponds to the electrical
   * current).
   */
  template<class Solver>
  class DOFVectorEvaluator {
    using ElementSearch = typename Solver::Traits::ElementSearch;
    using FunctionSpace = typename Solver::Traits::FunctionSpace;
    using GFS = typename FunctionSpace::GFS;
    using Scalar = typename Solver::Traits::CoordinateFieldType;
    using Matrix = DenseMatrix<Scalar>;
    using DomainDOFVector = Dune::PDELab::Backend::Vector<GFS, Scalar>;
    using GridView = typename Solver::Traits::GridView;
    enum { dim = GridView::dimension};
    static constexpr std::size_t OUTSIDE_DOMAIN = std::numeric_limits<std::size_t>::max();
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename Element::EntitySeed;
    using ElementGeometry = typename Element::Geometry;
    using GlobalCoordinate = typename ElementGeometry::GlobalCoordinate;
    using LocalCoordinate = typename ElementGeometry::LocalCoordinate;
  public:
    
    // Each DOF vector corresponds to a finite element trial function u.
    // "direct" means we evaluate u
    // "gradient" means we evaluate grad u
    // "current" means we evaluate - sigma * grad u, where sigma is the conductivity  
    // Note that we interprete the key "potential" to as meaning "direct".
    enum class EvaluationType {direct, gradient, current};
    
    // Some methods, such as the transfer matrix computation functions and the tDCS
    // matrix computation function, return the DOF vectors in a native form as actual vectors.
    // More concretely, they return matrices, where each row contains the coefficients defining a DOF vector.
    // This constructor is supposed to work with such a matrix.
    // Note that this class does not take ownership of the evaluationMatrix, and the user is responsible for 
    // ensuring that the object still exists when calling evaluate()
    explicit DOFVectorEvaluator(const Solver& solver, const Matrix& evaluationMatrix)
      : solver_(solver)
      , gfs_(solver_.functionSpace().getGFS())
      , elementSearch_(*(solver_.elementSearch()))
      , elementSeeds_(0)
      , localPositions_(0)
      , domainLabels_(0)
      , evaluationMatrixPtr_(&evaluationMatrix)
      , dofVectorPtrs_(0)
      , positionsBound_(false)
    {
      if(gfs_.ordering().size() != evaluationMatrix.cols()) {
        DUNE_THROW(Dune::Exception, "GridFunctionSpace has " << gfs_.ordering().size() << " degrees of freedom, but matrix has " << evaluationMatrix.cols() << " columns");
      }
    }
    
    // Note that this class does not take ownership of the DOFVector, and the user is responsible for 
    // ensuring that the object still exists when calling evaluate()
    explicit DOFVectorEvaluator(const Solver& solver, const DomainDOFVector& dofVector)
      : solver_(solver)
      , gfs_(solver.functionSpace().getGFS())
      , elementSearch_(*(solver_.elementSearch()))
      , elementSeeds_(0)
      , localPositions_(0)
      , domainLabels_(0)
      , evaluationMatrixPtr_(nullptr)
      , dofVectorPtrs_(1)
      , positionsBound_(false)
    {
      dofVectorPtrs_[0] = &dofVector;
    }
    
    // Note that this class does not take ownership of the DOFVectors, and the user is responsible for 
    // ensuring that the objects still exist when calling evaluate()
    explicit DOFVectorEvaluator(const Solver& solver, const std::vector<DomainDOFVector>& dofVectors)
      : solver_(solver)
      , gfs_(solver.functionSpace().getGFS())
      , elementSearch_(*(solver_.elementSearch()))
      , elementSeeds_(0)
      , localPositions_(0)
      , domainLabels_(0)
      , evaluationMatrixPtr_(nullptr)
      , dofVectorPtrs_(dofVectors.size())
      , positionsBound_(false)
    {
      for(std::size_t i = 0; i < dofVectors.size(); ++i) {
        dofVectorPtrs_[i] = &dofVectors[i];
      }
    }
    
    void bindPositions(const std::vector<GlobalCoordinate>& globalPositions)
    {
      std::size_t nrPositions = globalPositions.size();
      elementSeeds_.resize(nrPositions);
      localPositions_.resize(nrPositions);
      domainLabels_.resize(nrPositions);
      
      for(std::size_t i = 0; i < nrPositions; ++i) {
        auto search_result = elementSearch_.findEntity(globalPositions[i]);
        if(!search_result.has_value()) {
          DUNE_THROW(Dune::Exception, "coordinate is outside of the grid, or grid is not convex");
        }
        const auto& element = search_result.value();
        elementSeeds_[i] = element.seed();
        localPositions_[i] = element.geometry().local(globalPositions[i]);
        
        if constexpr(Solver::Traits::isFitted) {
          domainLabels_[i] = solver_.volumeConductor()->label(element);
        }
        else {
          std::optional<std::size_t> domainIndexOpt = solver_.domain().domainIndex(globalPositions[i]);
          if(domainIndexOpt.has_value()) {
            domainLabels_[i] = domainIndexOpt.value();
          }
          else {
            domainLabels_[i] = OUTSIDE_DOMAIN;
          }
        }
      }
      positionsBound_ = true;
      return;
    }

    std::unique_ptr<Matrix> evaluate(const Dune::ParameterTree& config)
    {
      if(!positionsBound_) {
        DUNE_THROW(Dune::Exception, "you need to bind positions before calling evaluate");
      }
      
      EvaluationType evaluationType = evaluationTypeFromString(config.get<std::string>("evaluation_return_type"));
      std::size_t nrValues = valuesPerPosition(evaluationType);
      
      std::size_t nrRows = evaluationMatrixPtr_ ? evaluationMatrixPtr_->rows() : dofVectorPtrs_.size();
      std::size_t nrPositions = elementSeeds_.size();
      
      std::unique_ptr<Matrix> output = std::make_unique<Matrix>(nrRows, nrValues * nrPositions);

#if HAVE_TBB
      int nr_threads = config.hasKey("numberOfThreads") ? config.get<int>("numberOfThreads") : tbb::task_arena::automatic;
      int grainSize = config.get<int>("grainSize", 1);
      tbb::task_arena arena(nr_threads);
      arena.execute([&]{
        tbb::parallel_for(
          tbb::blocked_range<std::size_t>(0, nrRows, grainSize),
          [&](const tbb::blocked_range<std::size_t>& range) {
            for(std::size_t i = range.begin(); i != range.end(); ++i) {
              if(evaluationMatrixPtr_) {
                DomainDOFVector wrappedVector(gfs_);
                extract_matrix_row(*evaluationMatrixPtr_, i, Dune::PDELab::Backend::native(wrappedVector));
                assembleOutputRow(wrappedVector, *output, i, evaluationType);
              }
              else {
                assembleOutputRow(*dofVectorPtrs_[i], *output, i, evaluationType);
              }
            }
          }
        );
      });
#else
      for(std::size_t i = 0; i < nrRows; ++i) {
        if(evaluationMatrixPtr_) {
          DomainDOFVector wrappedVector(gfs_);
          extract_matrix_row(*evaluationMatrixPtr_, i, Dune::PDELab::Backend::native(wrappedVector));
          assembleOutputRow(wrappedVector, *output, i, evaluationType);
        }
        else {
          assembleOutputRow(*dofVectorPtrs_[i], *output, i, evaluationType);
        }
      }
#endif
      return output;
    }

  private:
  
    // fitted case
    template<class SFINAEDummy = Solver>
    std::enable_if_t<Solver::Traits::isFitted && std::is_same_v<SFINAEDummy, Solver>> 
    assembleOutputRow(const DomainDOFVector& functionCoefficientVector, Matrix& output, std::size_t rowIndex, EvaluationType evaluationType) const
    {    
      if(evaluationType == EvaluationType::direct) {
        using DGF = Dune::PDELab::DiscreteGridFunction<GFS, DomainDOFVector>;
        using RangeType = typename DGF::Traits::RangeType;
        DGF function(gfs_, functionCoefficientVector);
        RangeType value;
        for(std::size_t i = 0; i < elementSeeds_.size(); ++i) {
          const auto& element = solver_.volumeConductor()->grid().entity(elementSeeds_[i]);
          function.evaluate(element, localPositions_[i], value);
          output(rowIndex, i) = value[0];
        }
      }
      else if(evaluationType == EvaluationType::gradient || evaluationType == EvaluationType::current){ 
        using DGFG = Dune::PDELab::DiscreteGridFunctionGradient<GFS, DomainDOFVector>;
        using GradientType = typename DGFG::Traits::RangeType;
        DGFG functionGradient(gfs_, functionCoefficientVector);
        GradientType gradient;
        for(std::size_t i = 0; i < elementSeeds_.size(); ++i) {
          const auto& element = solver_.volumeConductor()->grid().entity(elementSeeds_[i]);
          functionGradient.evaluate(element, localPositions_[i], gradient);
          if(evaluationType == EvaluationType::current) {
            GradientType currentDensity;
            solver_.volumeConductor()->tensor(element).mv(gradient, currentDensity);
            gradient = -currentDensity;
          }
          for(std::size_t j = 0; j < dim; ++j) {
            output(rowIndex, dim * i + j) = gradient[j];
          }
        } 
      }
      else {
        DUNE_THROW(Dune::NotImplemented, "evaluation type not implemented");
      }
      return;
    }
  
    // unfitted case
    template<class SFINAEDummy = Solver>
    std::enable_if_t<!Solver::Traits::isFitted && std::is_same_v<SFINAEDummy, Solver>>
    assembleOutputRow(const DomainDOFVector& functionCoefficientVector, Matrix& output, std::size_t rowIndex, EvaluationType evaluationType) const
    {      
      // In contrast to the fitted case, there does not seem to be an easy mechanism to wrap DOF vectors in GridFunctions, at least at the point of writing this function.
      // We thus directly assembly the output values from the finite element basis functions.
      using LocalFunctionSpace = typename Dune::PDELab::LocalFunctionSpace<GFS>;
      using LFSIndexCache = Dune::PDELab::LFSIndexCache<LocalFunctionSpace>;
      using ChildLFS = typename LocalFunctionSpace::template Child<0>::Type;
      using FESwitch = typename Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
      using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
      using GradientTypeMat = Dune::FieldMatrix<Scalar, 1, dim>;
      using GradientTypeVec = Dune::FieldVector<Scalar, dim>;
      using ValueTypeVec = Dune::FieldVector<Scalar, 1>;
      
      LocalFunctionSpace lfs(gfs_);
      LFSIndexCache lfsCache(lfs);
      
      if(evaluationType == EvaluationType::direct) {
        for(std::size_t i = 0; i < elementSeeds_.size(); ++i) {
          if(domainLabels_[i] == OUTSIDE_DOMAIN) {
            output(rowIndex, i) = 0.0;
            continue;
          }
          const auto& element = gfs_.gridView().grid().entity(elementSeeds_[i]);
          lfs.bind(element);
          lfsCache.update();
          std::vector<ValueTypeVec> basisValues;
          FESwitch::basis(lfs.child(domainLabels_[i]).finiteElement()).evaluateFunction(localPositions_[i], basisValues);
          output(rowIndex, i) = 0.0;
          for(std::size_t k = 0; k < lfs.child(domainLabels_[i]).size(); ++k) {
            output(rowIndex, i) += functionCoefficientVector[lfsCache.containerIndex(lfs.child(domainLabels_[i]).localIndex(k))] * basisValues[k][0];
          }
        }
      }
      else if(evaluationType == EvaluationType::gradient || evaluationType == EvaluationType::current) {
        // loop over all elements
        for(std::size_t i = 0; i < elementSeeds_.size(); ++i) {
          if(domainLabels_[i] == OUTSIDE_DOMAIN) {
            output(rowIndex, dim * i) = 0.0;
            output(rowIndex, dim * i + 1) = 0.0;
            output(rowIndex, dim * i + 2) = 0.0;
            continue;
          }
          const auto& element = gfs_.gridView().grid().entity(elementSeeds_[i]);
          lfs.bind(element);
          lfsCache.update();
          std::vector<GradientTypeMat> gradients;
          BasisSwitch::gradient(FESwitch::basis(lfs.child(domainLabels_[i]).finiteElement()), 
                                element.geometry(), 
                                localPositions_[i], 
                                gradients);
          GradientTypeVec gradient;
          gradient = 0.0;
          for(std::size_t k = 0; k < lfs.child(domainLabels_[i]).size(); ++k) {
            gradient += functionCoefficientVector[lfsCache.containerIndex(lfs.child(domainLabels_[i]).localIndex(k))] * gradients[k][0];
          }
          if(evaluationType == EvaluationType::current) {
            GradientTypeVec currentDensity;
            solver_.problem().A(domainLabels_[i]).mv(gradient, currentDensity);
            gradient = -currentDensity;
          }
          for(std::size_t j = 0; j < dim; ++j) {
            output(rowIndex, dim * i + j) = gradient[j];
          }
        }
      }
      else {
        DUNE_THROW(Dune::NotImplemented, "evaluation type not implemented");
      }
    }
  
    EvaluationType evaluationTypeFromString(const std::string& value)
    {
      if(value == "direct" || value == "potential") {
        return EvaluationType::direct;
      }
      else if(value == "gradient") {
        return EvaluationType::gradient;
      }
      else if(value == "current") {
        return EvaluationType::current;
      }
      else {
        DUNE_THROW(Dune::Exception, "invalid evaluation type (" << value << "), please choose one form {'direct', 'gradient', 'current'}");
      }
    }
    
    std::size_t valuesPerPosition(EvaluationType type)
    {
      std::size_t nrValues;
      if(type == EvaluationType::direct) {
        nrValues = 1;
      }
      else if(type == EvaluationType::gradient) {
        nrValues = dim;
      }
      else if(type == EvaluationType::current) {
        nrValues = dim;
      }
      else {
        DUNE_THROW(Dune::Exception, "unknown evaluation type");
      }
      
      return nrValues;
    }
  
    const Solver& solver_;
    const GFS& gfs_;
    const ElementSearch& elementSearch_;
    std::vector<ElementSeed> elementSeeds_;
    std::vector<LocalCoordinate> localPositions_;
    std::vector<std::size_t> domainLabels_;
    const Matrix* evaluationMatrixPtr_;
    std::vector<const DomainDOFVector*> dofVectorPtrs_;
    bool positionsBound_;
  };


} // namespace duneuro
#endif // DUNEURO_DOF_VECTOR_EVALUATOR_HH
