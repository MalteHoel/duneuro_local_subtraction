/*
 * uinfty.hh
 *
 * 	Implementation of the singularity potential u_infty and the gradient grad_u_infty as
 *analytic grid functions.
 *
 *  Created on: Jan 23, 2013
 *      Author: jakob
 */

#ifndef DUNEURO_SUBTRACTION_DG_UINFTY_HH
#define DUNEURO_SUBTRACTION_DG_UINFTY_HH

/**** includes ****/
#include <cmath>
#include <dune/pdelab/common/function.hh>
#include <dune/common/math.hh>

namespace duneuro
{
  template <typename GV, typename RF, class Enable = void>
  class InfinityPotential;

  /**** class definition ****/
  /*** u_infty ***/
  template <typename GV, typename RF>
  class InfinityPotential<GV, RF, typename std::enable_if<GV::dimension == 3>::type>
      : public Dune::PDELab::
            AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV, RF, 1>,
                                     InfinityPotential<GV, RF>>
  {
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, RF, 1> Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InfinityPotential<GV, RF>> B;

    /*** constructor(s) ***/
    explicit InfinityPotential(const GV& gv_) : B(gv_)
    {
      if (GV::dimension != 3) {
        DUNE_THROW(Dune::Exception, "expected dim == 3");
      }
    }

    /*** functions ***/
    /** evaluate the function for global coordinates **/
    inline void evaluateGlobal(const DomainType& x, RangeType& y) const
    {
      DomainType sigma_infty_inv_x_x_0;
      sigma_infty_inv.mv(x - x_0, sigma_infty_inv_x_x_0);
      RF numerator = M * sigma_infty_inv_x_x_0;
      RF denominator = sqrt((x - x_0) * sigma_infty_inv_x_x_0);
      denominator *= sqrt((x - x_0) * sigma_infty_inv_x_x_0);
      denominator *= sqrt((x - x_0) * sigma_infty_inv_x_x_0);

      RF ret = 1.0 / (4.0 * Dune::StandardMathematicalConstants<RF>::pi() * sqrt(sigma_infty.determinant()));
      ret *= numerator;
      ret /= denominator;

      y = ret;
    }

    /** set the parameters for the function **/
    void set_parameters(DomainType M_, DomainType x_0_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv_)
    {
      M = M_;
      x_0 = x_0_;
      sigma_infty = sigma_infty_;
      sigma_infty_inv = sigma_infty_inv_;
    }

  private:
    /*** dipole position, moment and sigma_infty ***/
    DomainType M;
    DomainType x_0;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv;
  };

  template <typename GV, typename RF, class Enable = void>
  class InfinityPotentialGradient;

  /*** grad_u_infty ***/
  template <typename GV, typename RF>
  class InfinityPotentialGradient<GV, RF, typename std::enable_if<GV::dimension == 3>::type>
      : public Dune::PDELab::
            AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV, RF,
                                                                              GV::dimension>,
                                     InfinityPotentialGradient<GV, RF>>
  {
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, RF, GV::dimension> Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InfinityPotentialGradient<GV, RF>> B;

    /*** constructor(s) ***/
    explicit InfinityPotentialGradient(const GV& gv_) : B(gv_)
    {
      if (GV::dimension != 3) {
        DUNE_THROW(Dune::Exception, "expected dim == 3");
      }
    }

    /** evaluate the gradient for global coordinates **/
    inline void evaluateGlobal(const DomainType& x, RangeType& y) const
    {
      RangeType siginvM;
      sigma_infty_inv.mv(M, siginvM);
      RangeType siginvx_x_0;
      sigma_infty_inv.mv(x - x_0, siginvx_x_0);
      y = siginvM;
      y *= (siginvx_x_0 * (x - x_0));
      RangeType y_help = siginvx_x_0;
      y_help *= 3 * (M * siginvx_x_0);
      y -= y_help;
      y /= 4 * Dune::StandardMathematicalConstants<RF>::pi() * sqrt(sigma_infty.determinant());
      y /= sqrt(siginvx_x_0 * (x - x_0));
      y /= siginvx_x_0 * (x - x_0);
      y /= siginvx_x_0 * (x - x_0);
    }

    /** set the parameters for the function **/
    void set_parameters(DomainType M_, DomainType x_0_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv_)
    {
      M = M_;
      x_0 = x_0_;
      sigma_infty = sigma_infty_;
      sigma_infty_inv = sigma_infty_inv_;
    }

  private:
    /*** dipole position, moment and sigma_infty ***/
    DomainType M;
    DomainType x_0;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv;
  };

  template <typename GV, typename RF>
  class InfinityPotential<GV, RF, typename std::enable_if<GV::dimension == 2>::type>
      : public Dune::PDELab::
            AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV, RF, 1>,
                                     InfinityPotential<GV, RF>>
  {
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, RF, 1> Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InfinityPotential<GV, RF>> B;

    /*** constructor(s) ***/
    explicit InfinityPotential(const GV& gv_) : B(gv_)
    {
      if (GV::dimension != 2) {
        DUNE_THROW(Dune::Exception, "expected dim == 2");
      }
    }

    /*** functions ***/
    /** evaluate the function for global coordinates **/
    inline void evaluateGlobal(const DomainType& x, RangeType& y) const
    {
      DomainType sigma_infty_inv_x_x_0;
      sigma_infty_inv.mv(x - x_0, sigma_infty_inv_x_x_0);
      RF numerator = M * sigma_infty_inv_x_x_0;
      RF denominator = (x - x_0) * sigma_infty_inv_x_x_0;

      RF ret = 1.0 / (2.0 * Dune::StandardMathematicalConstants<RF>::pi() * sqrt(sigma_infty.determinant()));
      ret *= numerator;
      ret /= denominator;
      y = ret;
    }

    /** set the parameters for the function **/
    void set_parameters(DomainType M_, DomainType x_0_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv_)
    {
      M = M_;
      x_0 = x_0_;
      sigma_infty = sigma_infty_;
      sigma_infty_inv = sigma_infty_inv_;
    }

  private:
    /*** dipole position, moment and sigma_infty ***/
    DomainType M;
    DomainType x_0;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv;
  };

  template <typename GV, typename RF>
  class InfinityPotentialGradient<GV, RF, typename std::enable_if<GV::dimension == 2>::type>
      : public Dune::PDELab::
            AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV, RF,
                                                                              GV::dimension>,
                                     InfinityPotentialGradient<GV, RF>>
  {
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, RF, GV::dimension> Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InfinityPotentialGradient<GV, RF>> B;

    /*** constructor(s) ***/
    explicit InfinityPotentialGradient(const GV& gv_) : B(gv_)
    {
      if (GV::dimension != 2) {
        DUNE_THROW(Dune::Exception, "expected dim == 2");
      }
    }

    /** evaluate the gradient for global coordinates **/
    inline void evaluateGlobal(const DomainType& x, RangeType& y) const
    {
      RangeType siginvM;
      sigma_infty_inv.mv(M, siginvM);
      RangeType siginvx_x_0;
      sigma_infty_inv.mv(x - x_0, siginvx_x_0);
      y = siginvM;
      RangeType y_help = siginvx_x_0;
      y_help *= 2. * (M * siginvx_x_0);
      y_help /= siginvx_x_0 * (x - x_0);
      y -= y_help;
      y /= 2. * Dune::StandardMathematicalConstants<RF>::pi() * sqrt(sigma_infty.determinant());
      y /= siginvx_x_0 * (x - x_0);
    }

    /** set the parameters for the function **/
    void set_parameters(DomainType M_, DomainType x_0_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_,
                        Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv_)
    {
      M = M_;
      x_0 = x_0_;
      sigma_infty = sigma_infty_;
      sigma_infty_inv = sigma_infty_inv_;
    }

  private:
    /*** dipole position, moment and sigma_infty ***/
    DomainType M;
    DomainType x_0;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty;
    Dune::FieldMatrix<RF, GV::dimension, GV::dimension> sigma_infty_inv;
  };
}

#endif
