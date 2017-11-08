#ifndef DUNEURO_FITTED_TENSOR_VTK_FUNCTOR_HH
#define DUNEURO_FITTED_TENSOR_VTK_FUNCTOR_HH

#include <dune/grid/io/file/vtk/function.hh>

#if HAVE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif

namespace duneuro
{
  /**
   * \brief vtk function representing the norm of a conductivity tensor
   *
   * Given a volume conductor, the vtk function evaluates the (piecewise constant) conductivity
   * tensor on each element and returns its infinity norm.
   */
  template <class VC>
  class FittedTensorNormFunctor : public Dune::VTKFunction<typename VC::GridView>
  {
  public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using Entity = typename GV::template Codim<0>::Entity;

    FittedTensorNormFunctor(std::shared_ptr<VC> volumeConductor) : volumeConductor_(volumeConductor)
    {
    }

    double evaluate(int, const Entity& e, const Dune::FieldVector<ctype, dim>&) const
    {
      return volumeConductor_->tensor(e).infinity_norm_real();
    }
    int ncomps() const
    {
      return 1;
    }
    std::string name() const
    {
      return "conductivity";
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
  };

#if HAVE_EIGEN
  /**
   * \brief vtk function representing a conductivity tensor
   *
   * Given a volume conductor, the vtk function evaluates the (piecewise constant) conductivity
   * tensor on each element. It computes the eigenvalues and corresponding eigenvectors. The
   * eigenvectors are scaled by their eigenvalues and one user-selected vector is returned as
   * the result of this functor.
   */
  template <class VC>
  class FittedTensorFunctor : public Dune::VTKFunction<typename VC::GridView>
  {
  public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using Entity = typename GV::template Codim<0>::Entity;

    FittedTensorFunctor(std::shared_ptr<VC> volumeConductor, unsigned int vectorIndex)
        : volumeConductor_(volumeConductor), cachedValues_(VC::dim), vectorIndex_(vectorIndex)
    {
    }

    double evaluate(int comp, const Entity& e, const Dune::FieldVector<ctype, dim>&) const
    {
      // cache scaled eigenvector so that it does not have to be recomputed between multiple calls
      // on the same element for different components
      if (cachedEntity_ != e) {
        // copy dune matrix to eigen format
        auto tensor = volumeConductor_->tensor(e);
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(VC::dim, VC::dim);
        for (unsigned int i = 0; i < VC::dim; ++i) {
          for (unsigned int j = 0; j < VC::dim; ++j) {
            matrix(i, j) = tensor[i][j];
          }
        }
        // scale eigen vectors and copy them to the output
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(matrix);
        Eigen::VectorXd v = es.eigenvectors().col(vectorIndex_);
        v *= es.eigenvalues()(vectorIndex_);
        for (unsigned int j = 0; j < VC::dim; ++j) {
          cachedValues_[j] = v(j);
        }
        cachedEntity_ = e;
      }
      return cachedValues_[comp];
    }

    int ncomps() const
    {
      return VC::dim;
    }

    std::string name() const
    {
      return std::string("conductivity_eigenvector_") + std::to_string(vectorIndex_);
    }

  private:
    std::shared_ptr<VC> volumeConductor_;
    mutable std::vector<double> cachedValues_;
    mutable Entity cachedEntity_;
    unsigned int vectorIndex_;
  };
#endif
}

#endif // DUNEURO_FITTED_TENSOR_VTK_FUNCTOR_HH
