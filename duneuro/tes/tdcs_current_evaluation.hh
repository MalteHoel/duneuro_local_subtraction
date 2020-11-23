#ifndef DUNEURO_TDCS_CURRENT_EVALUATION_HH
#define DUNEURO_TDCS_CURRENT_EVALUATION_HH
#include <duneuro/tes/tdcs_gradient_evaluation.hh>
namespace duneuro
{
template <typename GV, typename GFS>
class TDCSCurrentEvaluation : public TDCSGradientEvaluation<GV, GFS>{
    using Traits = TDCSGradientEvaluationTraits<GV,GFS>;
public:
  explicit TDCSCurrentEvaluation(const GFS& gfs) : TDCSGradientEvaluation<GV, GFS>(gfs) , gfs_(gfs) {}

  virtual std::vector<std::vector<double>> evaluateAtCoordinate(const typename Traits::Element& element,
     const typename Traits::LocalCoordinate& local, const DenseMatrix<double>& EvaluationMatrix, const typename Traits::SubTriangulation& subTriangulation) const override
    {
      typename Traits::ULFS ulfs(gfs_);
      typename Traits::UCache ucache(ulfs);
      typename Traits::UST ust(subTriangulation.gridView(), subTriangulation);
      std::vector<std::vector<double>> output(GV::dimension);   
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
          std::vector<Dune::FieldMatrix<typename Traits::Real, 1, GV::dimension>> gradphi(childLfs.size());
          Traits::FESwitch::basis(childLfs.finiteElement()).evaluateJacobian(local, gradphi);     // Gradient eval.
        //  const auto &conductivity = subTriangulation.tensor(element);        // still looking for easiest access to tensor
          for (std::size_t j = 0; j<GV::dimension; j++)
          {
            typename Traits::RangeDOFVector tmp(gfs_, 0.0);
            for (unsigned int i = 0; i < ucache.size(); ++i) {
                tmp[ucache.containerIndex(childLfs.localIndex(i))] += gradphi[i][0][j]; // gradphi is a Vector of 1xdim Matrices, so the first index
            }                                                                 // corresponds to the Ansatzfct. index, the third to the xyz-direction
         
          output[j] = matrix_dense_vector_product(EvaluationMatrix, Dune::PDELab::Backend::native(tmp));    // contains the derivative in j direction 
                                                                                                            // for each electrode
          for (auto& derivative : output[j])                              // multiply Gradient with cond. Tensor                                   
          {
              auto tmp2 = derivative;
              for (int k = 0; k< GV::dimension; k++)
              {
             //   derivative += tmp2*conductivity[j][k];
              }
          }                                                   
          }
          break;
          
      }
      return output;
    }
private:
const GFS& gfs_;

};
}
#endif // DUNEURO_TDCS_CURRENT_EVALUATION_HH