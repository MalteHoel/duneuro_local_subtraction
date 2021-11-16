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
        : gfs_(gfs), volumeConductor_(volumeConductor), current_(current), search_(gfs_.gridView())
    {
    }
    Dune::FieldVector<typename Traits::Real, GV::dimension>
    operator()(const typename Traits::Element& element,
               const typename Traits::LocalCoordinate localPos,
               const typename Traits::RangeDOFVector& coeffs) const
    {
      Dune::FieldVector<typename Traits::Real, GV::dimension> result;
      auto global = element.geometry().global(localPos);
      using DGFG = Dune::PDELab::DiscreteGridFunctionGradient<GFS, typename Traits::RangeDOFVector>;
      DGFG dgfg(gfs_, coeffs);

      // evaluate discrete grid function
      typename DGFG::Traits::RangeType y;
      dgfg.evaluate(element, localPos, y);
      if (current_) {
        const auto& conductivity = volumeConductor_.tensor(element);
        for (unsigned int i = 0; i < GV::dimension; ++i) {
          for (unsigned int j = 0; j < GV::dimension; ++j) {
            result[i] -= conductivity[i][j] * y[j];
          }
        }
      } else {
        for (unsigned int i = 0; i < GV::dimension; ++i) {
          result[i] -= y[i];
        }
      }
      return result;
    }

  private:
    const GFS& gfs_;
    const VC& volumeConductor_;
    const bool current_;
    KDTreeElementSearch<GV> search_;
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
    Dune::FieldVector<typename Traits::Real, GV::dimension>
    operator()(const typename Traits::Element& element,
               const typename Traits::LocalCoordinate localPos,
               const typename Traits::RangeDOFVector& coeffs) const
    {
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation_.gridView(), subTriangulation_);
      Dune::FieldVector<typename Traits::Real, GV::dimension> result;
      auto global = element.geometry().global(localPos);
      int comp = 0;

      Dune::FieldVector<typename Traits::Real, GV::dimension> y;
      std::vector<typename Traits::RangeType> phi; // storage for Ansatzfunction values

      ust.create(element); // Subtriangulate Element
      for (const auto& ep : ust) {
        const auto& geo = ep.geometry();
        const auto local = geo.local(global);
        const auto& refElement =
            Dune::ReferenceElements<typename GV::ctype, GV::dimension>::general(geo.type());
        if (!refElement.checkInside(local)) {
          continue;
        }
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
        Dune::FieldVector<typename Traits::Real, GV::dimension> y;
        y = 0;
        for (unsigned int i = 0; i < ucache.size(); ++i) {
          gradphi = 0;
          JgeoIT.umv(J[i][0], gradphi);
          // compute global gradient of shape function i (chain rule)
          const double& coeff = coeffs[ucache.containerIndex(childLfs.localIndex(i))];
          // sum up global gradients, weighting them with the appropriate coeff
          y.axpy(coeff, gradphi);
        }
        result = y;

        break;
      }
      return result;
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

  public:
    /**\brief constructor for the unfitted case
     */
    explicit TDCSGradientEvaluation(const GFS& gfs, const VC& volumeConductor, const bool current)
        : gradientEvaluator_(gfs, volumeConductor, current), current_(current), gfs_(gfs)
    {
    }

    virtual std::unique_ptr<DenseMatrix<double>>
    evaluate(const std::vector<typename Traits::Element>& elements,
             const std::vector<typename Traits::LocalCoordinate>& localPositions,
             const DenseMatrix<double>& EvaluationMatrix) const override
    {
      typename Traits::RangeDOFVector tmp(gfs_, 0.0);
      auto output = Dune::Std::make_unique<DenseMatrix<double>>(EvaluationMatrix.rows(),
                                                                3 * localPositions.size());
      for (unsigned int elec = 0; elec < EvaluationMatrix.rows(); ++elec) {
        tmp = 0.0;
        extract_matrix_row(EvaluationMatrix, elec, Dune::PDELab::Backend::native(tmp));
        std::size_t index = 0;
        std::vector<double> tmpVals(3 * localPositions.size());
        for (const auto& element : elements) {
          auto pot = gradientEvaluator_(elements[index], localPositions[index],
                                        tmp); // compute the gradient/current density
          for (unsigned int i = 0; i < GV::dimension; ++i) {
            tmpVals[3 * index + i] = (pot[i]);
          }
          index += 1;
        }
        set_matrix_row(*output, elec, tmpVals);
      }
      return output;
    }

  private:
    GradientEvaluator<GV, GFS, VC> gradientEvaluator_;
    const bool current_;
    const GFS& gfs_;
  };
}

#endif // DUNEURO_TDCS_GRADIENT_EVALUATION_HH