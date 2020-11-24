#ifndef DUNEURO_TDCS_POINT_EVALUATION_HH
#define DUNEURO_TDCS_POINT_EVALUATION_HH

#include <duneuro/tes/cutfem_tdcs_driver.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#include <duneuro/common/kdtree.hh>


namespace duneuro
{
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
  
  template <typename GV, typename GFS>
  class TDCSPointEvaluation : public TDCSEvaluationInterface<GV,GFS> {
    using Traits = TDCSPointEvaluationTraits<GV,GFS>;
  public:
  explicit TDCSPointEvaluation(const GFS& gfs) : gfs_(gfs) {}
#if HAVE_DUNE_UDG
  virtual std::vector<std::vector<double>> evaluate(const std::vector<typename Traits::Coordinate>& positions, const DenseMatrix<double>& EvaluationMatrix,
                                                    const typename Traits::SubTriangulation& subTriangulation) const override
    {
    std::vector<std::vector<double>> output(positions.size());
    KDTreeElementSearch<GV> search(gfs_.gridView());
    std::size_t index = 0;
    for (const auto& coord : positions) {
      
        const auto& element = search.findEntity(coord);
        if (!subTriangulation.isHostCell(element)) {
          DUNE_THROW(Dune::Exception, "element  at "
                                          << coord << " is not a host cell for any domain");
        }
        const auto localPos = element.geometry().local(coord);                      
        auto pot = evaluateAtCoordinate(element, localPos, EvaluationMatrix, subTriangulation); 
        output[index] = pot;
        index+=1;
    }
    return output;
    }
#endif
  private:
  const GFS& gfs_;

#if HAVE_DUNE_UDG
    std::vector<double> evaluateAtCoordinate(const typename Traits::Element& element,
     const typename Traits::LocalCoordinate& local, const DenseMatrix<double>& EvaluationMatrix, const typename Traits::SubTriangulation& subTriangulation) const
    {
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation.gridView(), subTriangulation);
      std::vector<double> output(EvaluationMatrix.rows());
      std::vector<typename Traits::RangeType> phi;                                 // storage for Ansatzfunction values     
      ust.create(element);                                        // splitting the Element 
      for (const auto& ep : ust) 
      {
          typename Traits::ChildLFS& childLfs(ulfs.child(ep.domainIndex() ) );     // chooses the correct Ansatzfunctionspace and binds it to the El.
          ulfs.bind(ep.subEntity(), true);
          ucache.update();

          if (childLfs.size() == 0)
          {
            continue;
          }
          Traits::FESwitch::basis(childLfs.finiteElement()).reset();

          phi.resize(childLfs.size());
          typename Traits::RangeDOFVector tmp(gfs_, 0.0);
          Traits::FESwitch::basis(childLfs.finiteElement()).evaluateFunction(local, phi);         // Ansatzfct eval.
          for (unsigned int i = 0; i < ucache.size(); ++i) {
              tmp[ucache.containerIndex(childLfs.localIndex(i))] += phi[i];
          }
          output = matrix_dense_vector_product(EvaluationMatrix, Dune::PDELab::Backend::native(tmp));
          break;
          
      }
      return output;
    }
#endif
  };
}

#endif // DUNEURO_TDCS_POINT_EVALUATION_HH