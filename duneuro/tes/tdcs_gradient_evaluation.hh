#ifndef DUNEURO_TDCS_GRADIENT_EVALUATION_HH
#define DUNEURO_TDCS_GRADIENT_EVALUATION_HH

/**
 *\file tdcs_gradient_evaluation.hh
 *\brief evaluates electric field at given positions, returns vector valued result as vector of whose 
 *       length is equal to the number of stimulation electrodes 
 */
#include <duneuro/tes/cutfem_tdcs_driver.hh>
#include <duneuro/tes/tdcs_evaluation_interface.hh>
#include <duneuro/common/kdtree.hh>
#include <duneuro/common/matrix_utilities.hh>

namespace duneuro
{
  #if HAVE_DUNE_UDG
template <typename GV, typename GFS, typename VC>
struct GradientEvaluator;
#endif
template <typename GV, typename GFS>
struct TDCSGradientEvaluationTraits{
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
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename ChildLFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
#endif
    using RangeType = typename BasisSwitch::Range;
    using RangeFieldType = typename BasisSwitch::RangeField;
    using Real = typename GV::ctype;
    using RangeDOFVector = Dune::PDELab::Backend::Vector<GFS, Real>;


};

  template <typename GV, typename GFS,typename VC>
  class TDCSGradientEvaluation : public TDCSEvaluationInterface<GV,GFS> {
      using Traits = TDCSGradientEvaluationTraits<GV,GFS>;

  public:
  /**\brief constructor for the unfitted case 
  */
    explicit TDCSGradientEvaluation(const GFS& gfs, const typename Traits::ST& subTriangulation, bool current)
                               : gradientEvaluator_ (gfs,subTriangulation), current_(current), gfs_(gfs){}

    virtual std::vector<std::vector<double>> evaluate(const std::vector<typename Traits::Element>& elements, const std::vector<typename Traits::LocalCoordinate>& localPositions,
                                                      const DenseMatrix<double>& EvaluationMatrix) const override
    {
      typename Traits::RangeDOFVector tmp(gfs_, 0.0);
      std::vector<std::vector<double>> output(localPositions.size());
      for (unsigned int elec = 0; elec<EvaluationMatrix.rows();++elec){
        tmp = 0.0;
        extract_matrix_row(EvaluationMatrix, elec, Dune::PDELab::Backend::native(tmp));
        std::size_t index = 0;
        for (const auto& element : elements) {                    
          auto pot = gradientEvaluator_(current_, elements[index], localPositions[index], tmp); // compute the potential
          for (unsigned int i = 0; i< GV::dimension;++i)
          {
            
            output[index].push_back(pot[i]); 
  
          }                                                 
          index+=1;
        }
      }
      return output;
    }
private:
GradientEvaluator<GV,GFS,VC> gradientEvaluator_;
bool current_;
const GFS& gfs_;
};
  #if HAVE_DUNE_UDG
  /**\class Functor that evaluates the elec. Field or Current density at a single point. Unfitted case.
 */
template <typename GV, typename GFS, typename ST>
struct GradientEvaluator
{
  using Traits = TDCSGradientEvaluationTraits<GV,GFS>;
  explicit GradientEvaluator (const GFS& gfs,const ST& subTriangulation) 
          : gfs_(gfs), subTriangulation_(subTriangulation) {}
  Dune::FieldVector<typename Traits::Real, GV::dimension> operator()(bool current, const typename Traits::Element& element, const typename Traits::LocalCoordinate localPos, const typename Traits::RangeDOFVector& coeffs) const
    {
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation_.gridView(), subTriangulation_);
      Dune::FieldVector<typename Traits::Real, GV::dimension> result;
      auto global = element.geometry().global(localPos); 
      int comp = 0;
      for (int i = 0;i<GV::dimension;i++)
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
      Dune::FieldVector<typename Traits::Real, GV::dimension> y;
      std::vector<typename Traits::RangeType> phi;                                 // storage for Ansatzfunction values

      ust.create(element);                                                          // Subtriangulate Element
      for (const auto& ep : ust) 
      {
        if(ep.domainIndex()!= comp)
        {continue;}
          typename Traits::ChildLFS& childLfs(ulfs.child(ep.domainIndex() ) );     // chooses the correct Ansatzfunctionspace and binds it to the El.
          ulfs.bind(ep.subEntity(), true);
          ucache.update();
          if (childLfs.size() == 0)
          {
            continue;
          }
          const auto JgeoIT = element.geometry().jacobianInverseTransposed(localPos);
        // get local Jacobians/gradients of the shape functions
        std::vector<Dune::FieldMatrix<typename Traits::Real, 1, GV::dimension>> J(childLfs.size());
        Traits::FESwitch::basis(childLfs.finiteElement()).reset();
        Traits::FESwitch::basis(childLfs.finiteElement()).evaluateJacobian(localPos, J);     // Gradients
        Dune::FieldVector<typename Traits::Real, GV::dimension> gradphi;
        Dune::FieldVector<typename Traits::Real, GV::dimension> y;

          y = 0; 
          for(unsigned int i = 0; i < ucache.size(); ++i) {
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
    const typename Traits::ST& subTriangulation_;
};
#endif
}


#endif // DUNEURO_TDCS_GRADIENT_EVALUATION_HH