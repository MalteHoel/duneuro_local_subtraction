#ifndef DUNEURO_Q1_GRADIENT_2D_HH
#define DUNEURO_Q1_GRADIENT_2D_HH

#include <dune/common/fmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace duneuro
{
  template <class D, class R>
  class Q1Gradient2DLocalBasis
  {
  public:
    typedef Dune::LocalBasisTraits<D, 2, Dune::FieldVector<D, 2>, R, 2, Dune::FieldVector<R, 2>,
                                   Dune::FieldMatrix<R, 2, 2>>
        Traits;

    //! \brief number of shape functions
    unsigned int size() const
    {
      return 3;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction(const typename Traits::DomainType& in,
                                 std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(3);
      out[0][0] = 1.0;
      out[0][1] = 0.0;
      out[1][0] = 0.0;
      out[1][1] = 1.0;
      out[2][0] = in[1];
      out[2][1] = in[0];
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian(const typename Traits::DomainType& in, // position
                     std::vector<typename Traits::JacobianType>& out) const // return value
    {
      out.resize(3);
      out[0][0][0] = 0;
      out[0][0][1] = 0;
      out[0][1][0] = 0;
      out[0][1][1] = 0;

      out[1][0][0] = 0;
      out[1][0][1] = 0;
      out[1][1][0] = 0;
      out[1][1][1] = 0;

      out[2][0][0] = 0;
      out[2][0][1] = 1;
      out[2][1][0] = 1;
      out[2][1][1] = 0;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial(const std::array<unsigned int, 2>& order,
                 const typename Traits::DomainType& in, // position
                 std::vector<typename Traits::RangeType>& out) const // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        auto const direction =
            std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        for (std::size_t i = 0; i < size(); ++i)
          out[i][0] = out[i][1] = 0;

        switch (direction) {
        case 0: out[2][1] = 1; break;
        case 1: out[2][0] = 1; break;
        default: DUNE_THROW(Dune::RangeError, "Component out of range.");
        }
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < 2; ++j)
            out[i][j] = 0;
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order() const
    {
      return 1;
    }
  };

  template <class LB, unsigned int size>
  class Q1Gradient2DLocalInterpolation
  {
  public:
    typedef typename LB::Traits::DomainType D;
    typedef typename LB::Traits::DomainFieldType DF;
    static const int dimD = LB::Traits::dimDomain;
    typedef typename LB::Traits::RangeType R;
    typedef typename LB::Traits::RangeFieldType RF;

    typedef Dune::QuadratureRule<DF, dimD> QR;

    Q1Gradient2DLocalInterpolation(const Dune::GeometryType& gt_, const LB& lb_)
        : gt(gt_), lb(lb_), Minv(0), qr(Dune::QuadratureRules<DF, dimD>::rule(gt, 2 * lb.order()))
    {
      // Compute inverse of the mass matrix of the local basis, and store it in Minv
      if (size != lb.size())
        DUNE_THROW(Dune::Exception, "size template parameter does not match size of local basis");

      for (const auto& qp : qr) {
        std::vector<R> base;
        lb.evaluateFunction(qp.position(), base);

        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            Minv[i][j] += qp.weight() * (base[i] * base[j]);
      }
      Minv.invert();
    }

    template <typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const
    {
      out.clear();
      out.resize(size, 0);

      for (const auto& qp : qr) {
        R y;
        f.evaluate(qp.position(), y);

        std::vector<R> base;
        lb.evaluateFunction(qp.position(), base);

        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            out[i] += Minv[i][j] * qp.weight() * (y * base[j]);
      }
    }

  private:
    Dune::GeometryType gt;
    const LB& lb;
    Dune::FieldMatrix<RF, size, size> Minv;
    const QR& qr;
  };

  class Q1Gradient2DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    Q1Gradient2DLocalCoefficients() : li(3)
    {
      for (std::size_t i = 0; i < 3; i++)
        li[i] = Dune::LocalKey(0, 0, i);
    }

    //! number of coefficients
    std::size_t size() const
    {
      return 3;
    }

    //! get i'th index
    const Dune::LocalKey& localKey(std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<Dune::LocalKey> li;
  };

  template <class D, class R>
  class Q1Gradient2DLocalFiniteElement
  {
    typedef Q1Gradient2DLocalBasis<D, R> LocalBasis;
    typedef Q1Gradient2DLocalCoefficients LocalCoefficients;
    typedef Q1Gradient2DLocalInterpolation<LocalBasis, 3> LocalInterpolation;

  public:
    /** \todo Please doc me !
     */
    typedef Dune::LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>
        Traits;

    /** \todo Please doc me !
     */
    Q1Gradient2DLocalFiniteElement()
      : gt(Dune::GeometryTypes::cube(2)), basis(), coefficients(), interpolation(gt, basis)
    {
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis() const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size() const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    Dune::GeometryType type() const
    {
      return gt;
    }

    Q1Gradient2DLocalFiniteElement* clone() const
    {
      return new Q1Gradient2DLocalFiniteElement(*this);
    }

  private:
    Dune::GeometryType gt;
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };
}

#endif // DUNEURO_Q1_GRADIENT_2D_HH
