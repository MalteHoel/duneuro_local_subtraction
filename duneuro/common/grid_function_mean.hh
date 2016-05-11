#ifndef DUNEURO_GRIDFUNCTIONMEAN_HH
#define DUNEURO_GRIDFUNCTIONMEAN_HH

#include <dune/common/ftraits.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/typetree/visitor.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/function/const.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class GFS>
  class GridFunctionMean : public Dune::TypeTree::TreeVisitor,
                           public Dune::TypeTree::DynamicTraversal
  {
  public:
    enum { dim = GFS::Traits::GridViewType::dimension };
    using RangeFieldType = double;
    using DOFVector = Dune::PDELab::Backend::Vector<GFS, RangeFieldType>;
    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    using IndexCache = Dune::PDELab::LFSIndexCache<LFS>;
    using FESwitch = Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;

    GridFunctionMean(std::shared_ptr<GFS> gfs, const Dune::ParameterTree& config)
        : gridFunctionSpace_(gfs), superIntegrationOrder_(config.get("superIntegrationOrder", 0))
    {
      assert(gridFunctionSpace_);
    }

    RangeFieldType evaluate(const DOFVector& x)
    {
      LFS lfs(*gridFunctionSpace_);
      IndexCache indexCache(lfs);

      RangeFieldType integral = 0.0;
      RangeFieldType volume = 0.0;

      std::vector<Dune::FieldVector<RangeFieldType, 1>> phi;

      for (const auto& element : elements(gridFunctionSpace_->gridView())) {
        lfs.bind(element);
        indexCache.update();
        phi.resize(indexCache.size());

        auto order = FESwitch::basis(lfs.finiteElement()).order();

        const auto& elementGeometry = element.geometry();

        using QRs = Dune::QuadratureRules<RangeFieldType, dim>;
        const auto& rule = QRs::rule(elementGeometry.type(), order);
        for (const auto& qp : rule) {
          FESwitch::basis(lfs.finiteElement()).evaluateFunction(qp.position(), phi);
          RangeFieldType u = 0.0;
          for (unsigned int i = 0; i < indexCache.size(); ++i) {
            u += x[indexCache.containerIndex(i)] * phi[i];
          }

          auto integrationElement = elementGeometry.integrationElement(qp.position());
          integral += qp.weight() * integrationElement * u;
          volume += qp.weight() * integrationElement;
        }
      }
      return integral / volume;
    }

  private:
    std::shared_ptr<GFS> gridFunctionSpace_;
    const unsigned int superIntegrationOrder_;
  };

  template <class GFS, class X>
  void subtract_mean_impl(const GFS& gfs, X& x,
                          const Dune::ParameterTree& config = Dune::ParameterTree(),
                          DataTree dataTree = DataTree())
  {
    Dune::Timer timer;
    GridFunctionMean<GFS> gfsMean(Dune::stackobject_to_shared_ptr(gfs), config);
    auto mean = gfsMean.evaluate(x);
    dataTree.set("mean_before_subtraction", mean);
    Dune::PDELab::ConstGridFunction<typename GFS::Traits::GridViewType, double> constGF(
        gfs.gridView(), mean);
    X mx(gfs, 0.0);
    Dune::PDELab::interpolate(constGF, gfs, mx);
    x -= mx;
    dataTree.set("mean_after_subtraction", gfsMean.evaluate(x));
    dataTree.set("time", timer.elapsed());
  }

  template <class Solver>
  void subtract_mean(const Solver& solver, typename Solver::Traits::DomainDOFVector& x,
                     const Dune::ParameterTree& config = Dune::ParameterTree(),
                     DataTree dataTree = DataTree())
  {
    subtract_mean_impl(solver.functionSpace().getGFS(), x, config, dataTree);
  }
}

#endif // DUNEURO_GRIDFUNCTIONMEAN_HH
