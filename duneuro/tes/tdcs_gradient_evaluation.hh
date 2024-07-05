#ifndef DUNEURO_TDCS_GRADIENT_EVALUATION_HH
#define DUNEURO_TDCS_GRADIENT_EVALUATION_HH

/**
 *\file tdcs_gradient_evaluation.hh
 *\brief evaluates electric field at given positions, returns vector valued result as vector of
 *whose length is equal to the number of stimulation electrodes
 */
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <duneuro/common/kdtree.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#if HAVE_DUNE_UDG
#include <dune/udg/pdelab/assembler/ulocalfunctionspace.hh>
#endif

namespace duneuro
{
#if HAVE_DUNE_UDG
  template <typename GV, typename GFS>
  struct UnfittedTDCSGradientEvaluationTraits {
    using BaseT = TDCSEvaluationInterface<GV, GFS>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using Real = typename GV::ctype;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
    using ChildLFS = typename ULFS::template Child<0>::Type;
    using FESwitch =
        Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using RangeType = typename BasisSwitch::Range;
  };
#endif
  template <typename GV, typename GFS>
  struct TDCSGradientEvaluationTraits {
    using BaseT = TDCSEvaluationInterface<GV, GFS>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using Real = typename GV::ctype;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;
  };

  /**\class Functor that evaluates the elec. Field or Current density at a single point. Fitted
   * case.
   */
  template <typename GV, typename GFS, typename VC>
  struct GradientEvaluator {
    using Traits = TDCSGradientEvaluationTraits<GV, GFS>;
    explicit GradientEvaluator(const GFS& gfs, const VC& volumeConductor, const bool current)
        : gfs_(gfs), volumeConductor_(volumeConductor), current_(current)
    {
    }
    
    virtual std::vector<Dune::FieldVector<typename Traits::Real, GV::dimension>>
    operator()(const typename Traits::Element& element,
               const typename Traits::LocalCoordinate localPos,
               const DenseMatrix<double>& EvaluationMatrix) const
    {
      std::vector<Dune::FieldVector<typename Traits::Real, GV::dimension>> result;
      result.resize(EvaluationMatrix.rows());
      auto global = element.geometry().global(localPos);
      typename Traits::RangeDOFVector coeffs(gfs_, 0.0);

      using DGFG = Dune::PDELab::DiscreteGridFunctionGradient<GFS, typename Traits::RangeDOFVector>;
      for (unsigned int nrElec = 0; nrElec < EvaluationMatrix.rows(); ++nrElec) {
        coeffs = 0.0;
        extract_matrix_row(EvaluationMatrix, nrElec, Dune::PDELab::Backend::native(coeffs));
      DGFG dgfg(gfs_, coeffs);
      // evaluate discrete grid function
      typename DGFG::Traits::RangeType y;
      dgfg.evaluate(element, localPos, y);
      if (current_) {
        const auto& conductivity = volumeConductor_.tensor(element);
        for (unsigned int i = 0; i < GV::dimension; ++i) {
          for (unsigned int j = 0; j < GV::dimension; ++j) {
            result[nrElec][i] -= conductivity[i][j] * y[j];
          }
        }
      } else {
        for (unsigned int i = 0; i < GV::dimension; ++i) {
          result[nrElec][i] -= y[i];
        }
      }
      }
      return result;
    }

  private:
    const GFS& gfs_;
    const VC& volumeConductor_;
    const bool current_;
  };

#if HAVE_DUNE_UDG
  /**\class Functor that evaluates the elec. Field or Current density at a single point. Unfitted
   * case.
   */
  template <typename GV, typename GFS>
  struct GradientEvaluator<GV, GFS, typename Dune::UDG::SimpleTpmcTriangulation<GV, GV>> {
    using Traits = UnfittedTDCSGradientEvaluationTraits<GV, GFS>;
    using ST = typename Dune::UDG::SimpleTpmcTriangulation<GV, GV>;
    explicit GradientEvaluator(const GFS& gfs, const ST& subTriangulation, const bool current)
        : gfs_(gfs), subTriangulation_(subTriangulation), current_(current)
    {
    }
    virtual std::vector<Dune::FieldVector<typename Traits::Real, GV::dimension>>
    operator()(const typename Traits::Element& element,
               const typename Traits::LocalCoordinate localPos,
               const DenseMatrix<double>& EvaluationMatrix) const
    {
      std::vector<Dune::FieldVector<typename Traits::Real, GV::dimension>> output;
      output.resize(EvaluationMatrix.rows());
      // check which domain the position belongs to by comparing it to the level set functions
      // first return to the global coordinate
      auto global = element.geometry().global(localPos);
      // access to the level set functions is given via the domainConfiguration
      const auto& dConf = subTriangulation_.domainConfiguration();
      // fill intRelPos with the info whether the point is in or outside the respective level set
      using InterfaceRelativePosition = Dune::UDG::InterfaceRelativePosition;
      std::vector<InterfaceRelativePosition> intRelPos;
      for (const auto& interface : subTriangulation_.domainConfiguration().interfaces() ){
        const auto& func = interface.function();
        auto val = func(global);
        if (val >0)
          intRelPos.push_back(InterfaceRelativePosition::exterior);
        else if (val<=0)
          intRelPos.push_back(InterfaceRelativePosition::interior);
        else
          intRelPos.push_back(InterfaceRelativePosition::any); // treat boundaries as "any"
      }
      // find the domain index that corresponds to intRelPos
      std::size_t index = 10000;
      for (const auto& domain : subTriangulation_.domainConfiguration().domains() ) {
        if (domain.contains(intRelPos)){
            index = domain.index();
            break;
        }
      }
      // if the point is outside all domains, return zeros
      if (index == 10000) {
        for (unsigned int nrElec = 0; nrElec< EvaluationMatrix.rows(); nrElec++){
         output[nrElec] = 0;  }
         return output;
      }
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation_.gridView(), subTriangulation_);
      auto result = Dune::FieldVector<typename Traits::Real, GV::dimension>(0) ;
      Dune::FieldVector<typename Traits::Real, GV::dimension> y;
      std::vector<typename Traits::RangeType> phi; // storage for Ansatzfunction values

      ust.create(element); // Subtriangulate Element
      for (const auto& ep : ust) {
        if (ep.domainIndex() != index)
          continue;
        typename Traits::ChildLFS& childLfs(ulfs.child(
            ep.domainIndex())); // chooses the correct Ansatzfunctionspace and binds it to the El.
        ulfs.bind(ep.subEntity(), true);
        ucache.update();
        if (childLfs.size() == 0) {
          continue;
        }
        const auto JgeoIT = element.geometry().jacobianInverseTransposed(localPos);
        // get local Jacobians/gradients of the shape functions
        std::vector<Dune::FieldMatrix<typename Traits::Real, 1, GV::dimension>> J(childLfs.size());
        Traits::FESwitch::basis(childLfs.finiteElement()).reset();
        Traits::FESwitch::basis(childLfs.finiteElement())
            .evaluateJacobian(localPos, J); // Gradients
        Dune::FieldVector<typename Traits::Real, GV::dimension> gradphi;
        //Dune::FieldVector<typename Traits::Real, GV::dimension> y;
        //y = 0;
        // EvaluationMatrix can no longer be accessed like a DOFVector using Multi-indices so we need flat indices
        using DOFVector = typename Traits::RangeDOFVector;
        using NativeVectorType = Dune::PDELab::Backend::Native<DOFVector>;
        // get the dof vector's block size 
        static constexpr int blockSize = NativeVectorType::block_type::dimension;
        for (unsigned int nrElec = 0; nrElec< EvaluationMatrix.rows(); nrElec++){
          output[nrElec] = 0;
          for (unsigned int i = 0; i < ucache.size(); ++i) {
            gradphi = 0;
            JgeoIT.umv(J[i][0], gradphi);
            // compute global gradient of shape function i (chain rule)
            auto CI = ucache.containerIndex(childLfs.localIndex(i));
            std::size_t coeffInd = CI[0]*blockSize + CI[1];
            const double& coeff = EvaluationMatrix(nrElec, coeffInd);
            // sum up global gradients, weighting them with the appropriate coeff

            output[nrElec].axpy(coeff, gradphi);
          }
        }
        return output;
      }
    }

  private:
    const GFS& gfs_;
    const ST& subTriangulation_;
    const bool current_;
  };
#endif

  template <typename GV, typename GFS, typename VC>
  class TDCSGradientEvaluation : public TDCSEvaluationInterface<GV, GFS>
  {
    using Traits = TDCSGradientEvaluationTraits<GV, GFS>;
    using RangeDOFVector = typename Traits::RangeDOFVector;
    using DGFG = Dune::PDELab::DiscreteGridFunctionGradient<GFS, RangeDOFVector>;
    using GradientType = typename DGFG::Traits::RangeType;

  public:
    /**\brief constructor for the fitted case
     */
    explicit TDCSGradientEvaluation(const GFS& gfs, const VC& volumeConductor, const bool current)
        : volumeConductor_(volumeConductor), current_(current), gfs_(gfs)
    {
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    evaluate(const std::vector<typename Traits::Element>& elements,
             const std::vector<typename Traits::LocalCoordinate>& localPositions,
             const DenseMatrix<double>& EvaluationMatrix, const Dune::ParameterTree& config) const override
    {
      std::size_t index = 0;
      auto output = std::make_unique<DenseMatrix<double>>(EvaluationMatrix.rows(),
                                                          3 * localPositions.size());

      // loop over all electrodes
#if HAVE_TBB
      int nr_threads = config.hasKey("numberOfThreads") ? config.get<int>("numberOfThreads") : tbb::task_arena::automatic;
      int grainSize = config.get<int>("grainSize", 1);
      tbb::task_arena arena(nr_threads);
      arena.execute([&]{
        tbb::parallel_for(
          tbb::blocked_range<std::size_t>(0, EvaluationMatrix.rows(), grainSize),
          [&](const tbb::blocked_range<std::size_t>& range) {
            for(size_t i = range.begin(); i != range.end(); ++i) {
              assembleGradientRow(elements, localPositions, EvaluationMatrix, i, *output);
            }
          }
        );
      });
#else
      for(size_t i = 0; i < EvaluationMatrix.rows(); ++i) {
        assembleGradientRow(elements, localPositions, EvaluationMatrix, i, *output);
      }
#endif

      return output;

    }

  private:
    
    void assembleGradientRow(const std::vector<typename Traits::Element>& elements,
                             const std::vector<typename Traits::LocalCoordinate>& localPositions,
                             const DenseMatrix<double>& EvaluationMatrix,
                             size_t rowIndex,
                             DenseMatrix<double>& output) const
    {      
      // get gradient current row
      RangeDOFVector tdcsPotential(gfs_);
      extract_matrix_row(EvaluationMatrix, rowIndex, Dune::PDELab::Backend::native(tdcsPotential));
      DGFG tdcsGradient(gfs_, tdcsPotential);
      GradientType gradient;
      
      // evaluate gradient at all positions specified in the arguments
      for(size_t j = 0; j < elements.size(); ++j) {
        tdcsGradient.evaluate(elements[j], localPositions[j], gradient);
        // correct to direction of current flow
        gradient *= -1.0;
        
        // potentially multiply with conductivity to get current
        if(current_) {
          GradientType currentDensity;
          volumeConductor_.tensor(elements[j]).mv(gradient, currentDensity);
          gradient = currentDensity;
        }
        
        // write gradient
        output(rowIndex, 3*j) = gradient[0];
        output(rowIndex, 3*j + 1) = gradient[1];
        output(rowIndex, 3*j + 2) = gradient[2];
      }
        
      return;  
    }
  
    const VC& volumeConductor_;
    const bool current_;
    const GFS& gfs_;
  };
  
#if HAVE_DUNE_UDG
template<typename GV, typename GFS>
  class TDCSGradientEvaluation<GV, GFS, Dune::UDG::SimpleTpmcTriangulation<GV, GV>>  : public TDCSEvaluationInterface<GV, GFS>
  {
    using Traits = TDCSGradientEvaluationTraits<GV, GFS>;

  public:
    /**\brief constructor for the fitted case
     */
    explicit TDCSGradientEvaluation(const GFS& gfs, const Dune::UDG::SimpleTpmcTriangulation<GV, GV>& volumeConductor, const bool current)
        : gradientEvaluator_(gfs, volumeConductor, current), current_(current), gfs_(gfs)
    {
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    evaluate(const std::vector<typename Traits::Element>& elements,
             const std::vector<typename Traits::LocalCoordinate>& localPositions,
             const DenseMatrix<double>& EvaluationMatrix, const Dune::ParameterTree& config) const override
    {
      std::size_t index = 0;
      auto output = std::make_unique<DenseMatrix<double>>(EvaluationMatrix.rows(),
                                                                3 * localPositions.size());
        for (const auto& element : elements) {
          auto pot = gradientEvaluator_(elements[index], localPositions[index],
                                        EvaluationMatrix); // compute the gradient/current density

          for (unsigned int i = 0; i<EvaluationMatrix.rows(); i++){
            (*output)(i,3*index)     = pot[i][0];
            (*output)(i,3*index + 1) = pot[i][1];
            (*output)(i,3*index + 2) = pot[i][2];
          }
          index++;
      }     

      return output;

    }

  private:
    GradientEvaluator<GV, GFS, Dune::UDG::SimpleTpmcTriangulation<GV, GV>> gradientEvaluator_;
    const bool current_;
    const GFS& gfs_;
  };
#endif
}

#endif // DUNEURO_TDCS_GRADIENT_EVALUATION_HH
