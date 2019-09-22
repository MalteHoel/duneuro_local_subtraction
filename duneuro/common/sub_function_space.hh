#ifndef DUNEURO_SUB_FUNCTION_SPACE_HH
#define DUNEURO_SUB_FUNCTION_SPACE_HH

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

namespace duneuro
{
  /*
   * \brief simple adapter class that adapts a function space on a sub entity set
   *
   * The input function space is assumed to be a scalar space.
   */
  template <class HostFS, class SubVC>
  class SubFunctionSpace
  {
  public:
    using ES = typename SubVC::EntitySet;
    using FEM = typename HostFS::FEM;
    using CONB = typename HostFS::CONB;
    using CON = typename HostFS::CON;
    using VBE = typename HostFS::VBE;
    using GFS = Dune::PDELab::GridFunctionSpace<ES, FEM, CON, VBE>;
    using DOF = Dune::PDELab::Backend::Vector<GFS, typename HostFS::NT>;
    using DGF = Dune::PDELab::DiscreteGridFunction<GFS, DOF>;
    using CC = typename GFS::template ConstraintsContainer<typename HostFS::NT>::Type;
    using NT = typename HostFS::NT;

    explicit SubFunctionSpace(const HostFS & hostfs, std::shared_ptr<SubVC> subVolumeConductor)
        : subVolumeConductor_(subVolumeConductor)
        , femp(std::make_shared<FEM>(hostfs.getFEM()))
        , gfsp(std::make_shared<GFS>(subVolumeConductor->entitySet(), femp))
        , ccp(std::make_shared<CC>())
    {
      gfsp->update();
    }

    FEM& getFEM()
    {
      return *femp;
    }
    const FEM& getFEM() const
    {
      return *femp;
    }

    GFS& getGFS()
    {
      return *gfsp;
    }

    const GFS& getGFS() const
    {
      return *gfsp;
    }
    CC& getCC()
    {
      return *ccp;
    }

    const CC& getCC() const
    {
      return *ccp;
    }

    template <class BCTYPE>
    void assembleConstraints(const BCTYPE& bctype, bool verbose = false)
    {
      ccp->clear();
      constraints(bctype, *gfsp, *ccp, verbose);
    }

    void clearConstraints()
    {
      ccp->clear();
    }

    void setConstrainedDOFS(DOF& x, typename HostFS::NT nt) const
    {
      set_constrained_dofs(*ccp, nt, x);
      // conb.make_consistent(*gfsp, x);
    }

    void setNonConstrainedDOFS(DOF& x, typename HostFS::NT nt) const
    {
      set_nonconstrained_dofs(*ccp, nt, x);
      // conb.make_consistent(*gfsp, x);
    }

    void copyConstrainedDOFS(const DOF& xin, DOF& xout) const
    {
      copy_constrained_dofs(*ccp, xin, xout);
      // conb.make_consistent(*gfsp, xout);
    }

    void copyNonConstrainedDOFS(const DOF& xin, DOF& xout) const
    {
      copy_nonconstrained_dofs(*ccp, xin, xout);
      // conb.make_consistent(*gfsp, xout);
    }

  private:
    std::shared_ptr<SubVC> subVolumeConductor_;
    // CONB conb;
    std::shared_ptr<FEM> femp;
    std::shared_ptr<GFS> gfsp;
    std::shared_ptr<CC> ccp;
  };
}

#endif // DUNEURO_SUB_FUNCTION_SPACE_HH
