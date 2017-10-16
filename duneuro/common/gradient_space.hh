#ifndef DUNEURO_GRADIENT_SPACE_HH
#define DUNEURO_GRADIENT_SPACE_HH

#include <dune/geometry/type.hh>

#include <dune/istl/solvercategory.hh>

#include <dune/localfunctions/monomial.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

#include <duneuro/common/p1gradientfem.hh>
#include <duneuro/common/q1gradientfem.hh>

namespace duneuro
{
  template <typename T, typename N, unsigned int degree, Dune::GeometryType::BasicType gt,
            Dune::SolverCategory::Category st = Dune::SolverCategory::sequential,
            typename VBET = Dune::PDELab::ISTL::
                VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,
                              T::dimension * Dune::QkStuff::QkSize<degree, T::dimension>::value>>
  class DGQkGradientSpace
  {
  public:
    // export types
    typedef T Grid;
    typedef typename T::LeafGridView GV;
    typedef typename T::ctype ctype;
    static const int dim = T::dimension;
    static const int dimworld = T::dimensionworld;
    typedef N NT;
    typedef duneuro::Q1GradientLocalFiniteElementMap<GV, NT, NT> FEM;
    typedef VBET VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints, VBE> GFS;
    typedef typename GFS::template ConstraintsContainer<N>::Type CC;
    using DOF = Dune::PDELab::Backend::Vector<GFS, N>;
    typedef Dune::PDELab::DiscreteGridFunction<GFS, DOF> DGF;
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> VTKF;

    // constructor making the grid function space an all that is needed
    DGQkGradientSpace(const GV& gridview) : gv(gridview), fem(), gfs(gv, fem), cc()
    {
      // initialize ordering
      gfs.update();
    }

    FEM& getFEM()
    {
      return fem;
    }
    const FEM& getFEM() const
    {
      return fem;
    }

    // return gfs reference
    GFS& getGFS()
    {
      return gfs;
    }

    // return gfs reference const version
    const GFS& getGFS() const
    {
      return gfs;
    }

    // return gfs reference
    CC& getCC()
    {
      return cc;
    }

    // return gfs reference const version
    const CC& getCC() const
    {
      return cc;
    }

  private:
    GV gv; // need this object here because FEM and GFS store a const reference !!
    FEM fem;
    GFS gfs;
    CC cc;
  };

  template <typename T, typename N, unsigned int degree, Dune::GeometryType::BasicType gt,
            Dune::SolverCategory::Category st = Dune::SolverCategory::sequential,
            typename VBET = Dune::PDELab::ISTL::
                VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,
                              T::dimension * Dune::MonomImp::Size<T::dimension, degree>::val>>
  class DGPkGradientSpace
  {
  public:
    // export types
    typedef T Grid;
    typedef typename T::LeafGridView GV;
    typedef typename T::ctype ctype;
    static const int dim = T::dimension;
    static const int dimworld = T::dimensionworld;
    typedef N NT;
    typedef duneuro::P1GradientLocalFiniteElementMap<GV, NT, NT> FEM;
    typedef VBET VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints, VBE> GFS;
    typedef typename GFS::template ConstraintsContainer<N>::Type CC;
    using DOF = Dune::PDELab::Backend::Vector<GFS, N>;
    typedef Dune::PDELab::DiscreteGridFunction<GFS, DOF> DGF;
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> VTKF;

    // constructor making the grid function space an all that is needed
    DGPkGradientSpace(const GV& gridview) : gv(gridview), fem(), gfs(gv, fem), cc()
    {
      // initialize ordering
      gfs.update();
    }

    FEM& getFEM()
    {
      return fem;
    }
    const FEM& getFEM() const
    {
      return fem;
    }

    // return gfs reference
    GFS& getGFS()
    {
      return gfs;
    }

    // return gfs reference const version
    const GFS& getGFS() const
    {
      return gfs;
    }

    // return gfs reference
    CC& getCC()
    {
      return cc;
    }

    // return gfs reference const version
    const CC& getCC() const
    {
      return cc;
    }

  private:
    GV gv; // need this object here because FEM and GFS store a const reference !!
    FEM fem;
    GFS gfs;
    CC cc;
  };
}
#endif // DUNEURO_GRADIENT_SPACE_HH
