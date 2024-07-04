#ifndef DUNEURO_TDCS_POTENTIAL_EVALUATION_HH
#define DUNEURO_TDCS_POTENTIAL_EVALUATION_HH
/**
 * \file tdcs_potential_evaluation.hh
 * \brief evaluates electric potential at given points
 */
#include <duneuro/common/kdtree.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>

namespace duneuro
{
#if HAVE_DUNE_UDG

#if HAVE_DUNE_UDG
  template <typename GV, typename GFS, typename VC>
  struct PotentialEvaluator;
#endif
  template <typename GV, typename GFS>
  struct TDCSPotentialEvaluationTraits {
    using BaseT = TDCSEvaluationInterface<GV, GFS>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using Coordinate = typename BaseT::Coordinate;
    using ST = typename BaseT::SubTriangulation;
#if HAVE_DUNE_UDG
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
    using ChildLFS = typename ULFS::template Child<0>::Type;
    using FESwitch =
        Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
#endif
    using RangeType = typename BasisSwitch::Range;
    using RangeFieldType = typename BasisSwitch::RangeField;
    using Real = typename GV::ctype;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;
  };
  /**\class Implementation of the EvaluationInterface for Potential Evaluation at given Points
   * @tparam VC should be either the SubTriangulation for Unfitted or a VolumeConductor for fitted
   */
  template <typename GV, typename GFS, typename VC>
  class TDCSPotentialEvaluation : public TDCSEvaluationInterface<GV, GFS>
  {
    using Traits = TDCSPotentialEvaluationTraits<GV, GFS>;

  public:
    /**\brief constructor for the unfitted case
     */
    explicit TDCSPotentialEvaluation(const GFS& gfs, const typename Traits::ST& subTriangulation)
        : potentialEvaluator_(gfs, subTriangulation)
    {
    }

    virtual std::vector<std::vector<double>>
    evaluate(const std::vector<typename Traits::Element>& elements,
             const std::vector<typename Traits::LocalCoordinate>& localPositions,
             const DenseMatrix<double>& EvaluationMatrix) const override
    {
      std::vector<std::vector<double>> output(localPositions.size());

      std::size_t index = 0;
      for (const auto& element : elements) {
        auto pot = potentialEvaluator_(elements[index], localPositions[index],
                                       EvaluationMatrix); // compute the potential
        output[index] = pot;
        index += 1;
      }
      return output;
    }

  private:
    PotentialEvaluator<GV, GFS, VC> potentialEvaluator_;
  };
#if HAVE_DUNE_UDG
  /**\class Functor that evaluates the elec. Pot. at a single point. Unfitted case.
   */
  template <typename GV, typename GFS, typename ST>
  struct PotentialEvaluator {
    using Traits = TDCSPotentialEvaluationTraits<GV, GFS>;
    explicit PotentialEvaluator(const GFS& gfs, const ST& subTriangulation)
        : gfs_(gfs), subTriangulation_(subTriangulation)
    {
    }
    std::vector<double> operator()(const typename Traits::Element& element,
                                   const typename Traits::LocalCoordinate localPos,
                                   const DenseMatrix<double>& EvaluationMatrix) const
    {
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation_.gridView(), subTriangulation_);
      std::vector<double> output(EvaluationMatrix.rows());
      std::vector<typename Traits::RangeType> phi; // storage for Ansatzfunction values
      ust.create(element); // splitting of the Element
      for (const auto& ep : ust) {
        typename Traits::ChildLFS& childLfs(ulfs.child(
            ep.domainIndex())); // chooses the correct Ansatzfunctionspace and binds it to the El.
        ulfs.bind(ep.subEntity(), true);
        ucache.update();

        if (childLfs.size() == 0) {
          continue;
        }
        Traits::FESwitch::basis(childLfs.finiteElement())
            .reset(); // reset from previous computations

        phi.resize(childLfs.size());
        typename Traits::RangeDOFVector tmp(gfs_, 0.0);
        Traits::FESwitch::basis(childLfs.finiteElement())
            .evaluateFunction(localPos, phi); // Ansatzfct eval.
        for (unsigned int i = 0; i < ucache.size(); ++i) {
          tmp[ucache.containerIndex(childLfs.localIndex(i))] += phi[i];
        }
        // tmp contains the Ansatzfct. values at the given position, multiplication with coefficient
        // matrix yields vector containing the respective potential value for each stimulation
        // electrode
        output = matrix_dense_vector_product(EvaluationMatrix, Dune::PDELab::Backend::native(tmp));
        break;
      }
      return output;
    }

  private:
    const GFS& gfs_;
    const typename Traits::ST& subTriangulation_;
  };
#endif

#endif // HAVE_DUNE_UDG
}

#endif // DUNEURO_TDCS_POTENTIAL_EVALUATION_HH
