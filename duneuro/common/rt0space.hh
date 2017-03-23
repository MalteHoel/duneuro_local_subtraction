#ifndef DUNEURO_RT0SPACE_HH
#define DUNEURO_RT0SPACE_HH

#include <dune/geometry/type.hh>

#include <dune/istl/solvercategory.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

namespace duneuro
{
  template <typename T, typename N, unsigned int degree, Dune::GeometryType::BasicType gt,
            Dune::SolverCategory::Category st = Dune::SolverCategory::sequential,
            typename VBET = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>>
  class RT0Space
  {
  public:
    // export types
    typedef T Grid;
    typedef typename T::LeafGridView GV;
    typedef typename T::ctype ctype;
    static const int dim = T::dimension;
    static const int dimworld = T::dimensionworld;
    typedef N NT;
    typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV, N, N, 0, gt> FEM;
    typedef VBET VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints, VBE> GFS;
    typedef typename GFS::template ConstraintsContainer<N>::Type CC;
    using DOF = Dune::PDELab::Backend::Vector<GFS, N>;
    typedef Dune::PDELab::DiscreteGridFunction<GFS, DOF> DGF;
    typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> VTKF;

    // constructor making the grid function space an all that is needed
    RT0Space(const GV& gridview) : gv(gridview)
    {
      femp = std::shared_ptr<FEM>(new FEM(gv));
      gfsp = std::shared_ptr<GFS>(new GFS(gv, *femp));
      // initialize ordering
      gfsp->update();
      ccp = std::shared_ptr<CC>(new CC());
    }

    FEM& getFEM()
    {
      return *femp;
    }
    const FEM& getFEM() const
    {
      return *femp;
    }

    // return gfs reference
    GFS& getGFS()
    {
      return *gfsp;
    }

    // return gfs reference const version
    const GFS& getGFS() const
    {
      return *gfsp;
    }

    // return gfs reference
    CC& getCC()
    {
      return *ccp;
    }

    // return gfs reference const version
    const CC& getCC() const
    {
      return *ccp;
    }

  private:
    GV gv; // need this object here because FEM and GFS store a const reference !!
    std::shared_ptr<FEM> femp;
    std::shared_ptr<GFS> gfsp;
    std::shared_ptr<CC> ccp;
  };
}
#endif // DUNEURO_RT0SPACE_HH
