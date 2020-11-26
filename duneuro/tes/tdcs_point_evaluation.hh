#ifndef DUNEURO_TDCS_POINT_EVALUATION_HH
#define DUNEURO_TDCS_POINT_EVALUATION_HH

#include <duneuro/tes/cutfem_tdcs_driver.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#include <duneuro/common/kdtree.hh>

// evaluates electric potential at given points
namespace duneuro
{
#if HAVE_DUNE_UDG
template <typename GV, typename GFS, typename VC>
struct TDCSCoordinatePointEvaluation;
#endif
template<typename GV, typename GFS>
struct TDCSPointEvaluationTraits{
    using BaseT = TDCSEvaluationInterface<GV, GFS>;
    using LocalCoordinate = typename BaseT::LocalCoordinate;
    using Element = typename BaseT::Element;
    using Coordinate = typename BaseT::Coordinate;
    using SubTriangulation = typename BaseT::SubTriangulation;
#if HAVE_DUNE_UDG
    using ULFS = Dune::PDELab::UnfittedLocalFunctionSpace<GFS>;
    using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;
    using UCache = Dune::PDELab::LFSIndexCache<ULFS>;
    using ChildLFS = typename ULFS::template Child<0>::Type;
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
#endif
    using RangeType = typename BasisSwitch::Range;
    using RangeFieldType = typename BasisSwitch::RangeField;
    using Real = typename GV::ctype;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;
};

template <typename GV, typename GFS, typename VC>
class TDCSPointEvaluation : public TDCSEvaluationInterface<GV,GFS> {
    using Traits = TDCSPointEvaluationTraits<GV,GFS>;
  public:
    explicit TDCSPointEvaluation(const GFS& gfs, const typename Traits::SubTriangulation& subTriangulation)
                               : evaluatorBackend_ (gfs,subTriangulation) {}

    virtual std::vector<std::vector<double>> evaluate(const std::vector<typename Traits::Coordinate>& positions,
                                                      const DenseMatrix<double>& EvaluationMatrix) const override
    {

      std::vector<std::vector<double>> output(positions.size());

      std::size_t index = 0;
      for (const auto& coord : positions) {                    
        auto pot = evaluatorBackend_(coord, EvaluationMatrix); // compute the potential
        output[index] = pot;                                                      
        index+=1;
      }
      return output;
    }
private:
TDCSCoordinatePointEvaluation<GV,GFS,VC> evaluatorBackend_;
};
#if HAVE_DUNE_UDG
template <typename GV, typename GFS, typename SubTriangulation>
struct TDCSCoordinatePointEvaluation
{
  using Traits = TDCSPointEvaluationTraits<GV,GFS>;
  explicit TDCSCoordinatePointEvaluation (const GFS& gfs,const SubTriangulation& subTriangulation) 
          : gfs_(gfs), subTriangulation_(subTriangulation) {}
  std::vector<double> operator()(const typename Traits::BaseT::Coordinate coord, const DenseMatrix<double>& EvaluationMatrix) const
    {
        KDTreeElementSearch<GV> search(gfs_.gridView());
        const auto& element = search.findEntity(coord);       // find cut-cell containing the position
        if (!subTriangulation_.isHostCell(element)) {
          DUNE_THROW(Dune::Exception, "element  at "
                                          << coord << " is not a host cell for any domain");
        }
        const auto localPos = element.geometry().local(coord);  
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation_.gridView(), subTriangulation_);
      std::vector<double> output(EvaluationMatrix.rows());
      std::vector<typename Traits::RangeType> phi;                                  // storage for Ansatzfunction values     
      ust.create(element);                                                          // splitting of the Element 
      for (const auto& ep : ust) 
      {
          typename Traits::ChildLFS& childLfs(ulfs.child(ep.domainIndex() ) );     // chooses the correct Ansatzfunctionspace and binds it to the El.
          ulfs.bind(ep.subEntity(), true);
          ucache.update();

          if (childLfs.size() == 0)
          {
            continue;
          }
          Traits::FESwitch::basis(childLfs.finiteElement()).reset();                  // reset from previous computations

          phi.resize(childLfs.size());
          typename Traits::RangeDOFVector tmp(gfs_, 0.0);
          Traits::FESwitch::basis(childLfs.finiteElement()).evaluateFunction(localPos, phi);         // Ansatzfct eval.
          for (unsigned int i = 0; i < ucache.size(); ++i) {
              tmp[ucache.containerIndex(childLfs.localIndex(i))] += phi[i];                   
          }
      // tmp contains the Ansatzfct. values at the given position, multiplication with coefficient matrix yields 
      // vector containing the respective potential value for each stimulation electrode
          output = matrix_dense_vector_product(EvaluationMatrix, Dune::PDELab::Backend::native(tmp));
          break;
          
      }
      return output;
    }
  private:
    const GFS& gfs_;
    const typename Traits::SubTriangulation& subTriangulation_;
  };
#endif
}

#endif // DUNEURO_TDCS_POINT_EVALUATION_HH