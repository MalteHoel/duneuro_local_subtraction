#ifndef DUNEURO_FITTED_TENSOR_VTK_FUNCTOR_HH
#define DUNEURO_FITTED_TENSOR_VTK_FUNCTOR_HH

#include <cmath>
#include <algorithm>

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

    FittedTensorNormFunctor(std::shared_ptr<const VC> volumeConductor)
        : volumeConductor_(volumeConductor)
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
    std::shared_ptr<const VC> volumeConductor_;
  };

  /**
   * \brief vtk function representing a complete conductivity tensor
   *
   * Given a fitted volume conductor, this function evaluates the (piecewise constant) conductivity 
   * tensor on each element and then forwards it to the VTK writer.
   */
   template <class VC>
   class FittedTensorFunctor : public Dune::VTKFunction<typename VC::GridView>
   {
   public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum {dim = GV::dimension};
    using Entity = typename GV::template Codim<0>::Entity;
    
    FittedTensorFunctor(std::shared_ptr<const VC> volumeConductor)
      : volumeConductor_(volumeConductor)
    {
    }
    
    double evaluate(int component, const Entity& e, const Dune::FieldVector<ctype, dim>&) const
    {
      int row = component / dim;
      int col = component % dim;
      auto tensor = volumeConductor_->tensor(e);
      return tensor[row][col];
    }
    
    int ncomps() const
    {
      return dim * dim;
    }
    
    std::string name() const
    {
      return "conductivity_tensor";
    }
   
   private:
    std::shared_ptr<const VC> volumeConductor_;
   };
   
 /**
   * \brief vtk function representing the fractional anisotropy
   *
   * Given a fitted volume conductor, this function evaluates the (piecewise constant) conductivity 
   * tensor on each element and then computes its fractional anisotropy.
   */ 
   template <class VC>
   class FittedTensorFractionalAnisotropyFunctor : public Dune::VTKFunction<typename VC::GridView>
   {
   public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum {dim = GV::dimension};
    using Entity = typename GV::template Codim<0>::Entity;
    
    FittedTensorFractionalAnisotropyFunctor(std::shared_ptr<const VC> volumeConductor)
      : volumeConductor_(volumeConductor)
    {
    }
    
    /* Let D be a nonzero symmetric positive semidefinite tensor, with eigenvalues lambda_1, lambda_2, and lambda_3. Then The fractional anisotropy can be computed via
     *  FA(D) = sqrt((1/2) * (3 - Trace(D)^2/Trace(D^2)))
     * see e.g. https://en.wikipedia.org/wiki/Fractional_anisotropy.
     * Then, 0 <= FA(D) <= 1, with 0 = FA(D) if, and only if, D is a multiple of the identity matrix (i.e. isotropic), and FA(D) = 1 if, and only if, exactly one of the eigenvalues
     * is non-zero.
     */
    double evaluate(int, const Entity& e, const Dune::FieldVector<ctype, dim>&) const
    {
      auto tensor = volumeConductor_->tensor(e);
      auto squared_tensor = tensor * tensor;
      
      double tensor_trace = 0.0;
      double squared_tensor_trace = 0.0;
      
      for(size_t i = 0; i < dim; ++i) {
        tensor_trace += tensor[i][i];
        squared_tensor_trace += squared_tensor[i][i]; 
      }
      
      double fa_squared = 0.5 * (3.0 - (tensor_trace * tensor_trace)/squared_tensor_trace);
      
      // in exact arithmetic, one has 0 <= FA <= 1. Due to rounding errors, we might go slightly
      // outside of this range. To prevent this, we clip the squared FA to the interval [0, 1].
      double clipped_fa_squared = std::max(std::min(fa_squared, 1.0), 0.0);
      
      return std::sqrt(clipped_fa_squared);
    }
    
    int ncomps() const
    {
      return 1;
    }
    
    std::string name() const
    {
      return "fractional_anisotropy";
    }
   
   private:
    std::shared_ptr<const VC> volumeConductor_;
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
  class FittedTensorEigenvectorFunctor : public Dune::VTKFunction<typename VC::GridView>
  {
  public:
    using GV = typename VC::GridView;
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using Entity = typename GV::template Codim<0>::Entity;

    FittedTensorEigenvectorFunctor(std::shared_ptr<const VC> volumeConductor, unsigned int vectorIndex)
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
    std::shared_ptr<const VC> volumeConductor_;
    mutable std::vector<double> cachedValues_;
    mutable Entity cachedEntity_;
    unsigned int vectorIndex_;
  };
#endif
}

#endif // DUNEURO_FITTED_TENSOR_VTK_FUNCTOR_HH
