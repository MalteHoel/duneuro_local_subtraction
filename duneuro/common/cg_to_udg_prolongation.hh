#ifndef DUNEURO_CG_TO_UDG_PROLONGATION_HH
#define DUNEURO_CG_TO_UDG_PROLONGATION_HH

#include <memory>

#include <dune/common/fmatrix.hh>

#include <dune/grid/common/scsgmapper.hh>

#include <dune/istl/bcrsmatrix.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>

namespace duneuro
{
  namespace cg2udgdetail
  {
    template <class FE>
    struct ShapeFunctionOnBoundingBoxTraits {
      typedef Dune::FiniteElementInterfaceSwitch<FE> FESwitch;
      typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
      typedef typename BasisSwitch::DomainField DF;
      typedef typename BasisSwitch::Range RangeType;
      enum { dim = BasisSwitch::dimDomainLocal };

    };
    // evaluate a localfunction in local coordinates of a bounding box
    template <class FE, class BB, class E>
    class ShapeFunctionOnBoundingBox
    {

      const FE& fe_;
      const int comp_;
      const BB& boundingBox_;
      const E& element_;

    public:
      using Traits = ShapeFunctionOnBoundingBoxTraits<FE>;

      ShapeFunctionOnBoundingBox(const FE& fe, int comp, const BB& boundingBox, const E& element)
          : fe_(fe), comp_(comp), boundingBox_(boundingBox), element_(element)
      {
      }

      void evaluate(const Dune::FieldVector<typename Traits::DF, Traits::dim>& x,
                    typename Traits::RangeType& y) const
      {
        std::vector<typename Traits::RangeType> v;
        fe_.localBasis().evaluateFunction(element_.geometry().local(boundingBox_.global(x)), v);
        y = v[comp_];
      }
    };

    template <class FE, class BB, class E>
    std::unique_ptr<ShapeFunctionOnBoundingBox<FE, BB, E>>
    make_shape_function_on_bounding_box(const FE& fe, int comp, const BB& boundingBox,
                                        const E& element)
    {
      return std::make_unique<ShapeFunctionOnBoundingBox<FE, BB, E>>(fe, comp, boundingBox,
                                                                           element);
    }

    template <class E, class IC, class CI, class M>
    class ComputeCG2UDGVisitor : public Dune::TypeTree::DefaultVisitor,
                                 public Dune::TypeTree::DynamicTraversal,
                                 public Dune::TypeTree::VisitTree
    {
      M& mat_;
      const E& element_;
      const IC& indexCache_;
      const CI& cutCellInformation_;
      std::vector<std::size_t> vertexIndices_;
      const std::map<std::pair<std::size_t, std::size_t>, std::size_t>& vertexAndDomainToLinear_;

    public:
      template <class GV>
      ComputeCG2UDGVisitor(
          M& mat, const E& element, const IC& indexCache, const CI& cutCellInformation,
          const GV& gridView,
          const std::map<std::pair<std::size_t, std::size_t>, std::size_t>& vertexAndDomainToLinear)
          : mat_(mat)
          , element_(element)
          , indexCache_(indexCache)
          , cutCellInformation_(cutCellInformation)
          , vertexAndDomainToLinear_(vertexAndDomainToLinear)
      {
        Dune::SingleCodimSingleGeomTypeMapper<GV, GV::dimension> mapper(gridView);
        for (unsigned int i = 0; i < element.geometry().corners(); ++i) {
          vertexIndices_.push_back(mapper.subIndex(element, i, GV::dimension));
        }
      }

      template <class LFS, class TreePath>
      void leaf(const LFS& lfs, TreePath treePath) const
      {
        if (lfs.size() == 0)
          return;
        typedef typename LFS::Traits::FiniteElementType UDGFEM;
        typedef Dune::FiniteElementInterfaceSwitch<UDGFEM> FESwitch;
        typedef Dune::BasisInterfaceSwitch<typename FESwitch::Basis> BasisSwitch;
        typedef typename BasisSwitch::DomainField DF;
        Dune::PDELab::QkDGLocalFiniteElementMap<double, double, 1, LFS::Traits::GridFunctionSpace::
                                                                       Traits::GridView::dimension>
            cgfem;
        const auto& fe = cgfem.find(element_);
        const auto& bbox =
            cutCellInformation_.information(element_, treePath.element(0)).boundingBox;
        std::vector<DF> v;
        for (unsigned int i = 0; i < fe.localBasis().size(); ++i) {
          auto sf = make_shape_function_on_bounding_box(fe, i, bbox, element_);
          lfs.finiteElement().localFiniteElement().localInterpolation().interpolate(*sf, v);
          auto it = vertexAndDomainToLinear_.find({vertexIndices_[i], treePath.element(0)});
          if (it == vertexAndDomainToLinear_.end()) {
            DUNE_THROW(Dune::Exception, "vertex " << vertexIndices_[i] << " domain "
                                                  << treePath.element(0) << " not found in map");
          }
          for (unsigned int j = 0; j < v.size(); ++j) {
            const auto& index = indexCache_.containerIndex(lfs.localIndex(j));
            if (!mat_.exists(index[1], it->second)) {
              DUNE_THROW(Dune::Exception, "index (" << index[1] << "," << it->second
                                                    << ") out of bounds");
            }
            mat_[index[1]][it->second][j][0] = v[j];
          }
        }
      }
    };

    template <class E, class IC, class CI, class M, class GV>
    std::unique_ptr<ComputeCG2UDGVisitor<E, IC, CI, M>> make_compute_cg_2_udg_visitor(
        M& mat, const E& element, const IC& indexCache, const CI& cutCellInformation,
        const GV& gridView,
        const std::map<std::pair<std::size_t, std::size_t>, std::size_t>& vertexAndDomainToLinear)
    {
      return std::make_unique<ComputeCG2UDGVisitor<E, IC, CI, M>>(
          mat, element, indexCache, cutCellInformation, gridView, vertexAndDomainToLinear);
    }
  }

  template <class GFS, class ST>
  std::unique_ptr<Dune::BCRSMatrix<Dune::FieldMatrix<double, Dune::QkStuff::
                                                                 QkSize<1, GFS::Traits::GridView::
                                                                               dimension>::value,
                                                     1>>>
  compute_cg_to_udg_prolongation(const GFS& gfs, const ST& st)
  {
    Dune::SingleCodimSingleGeomTypeMapper<typename GFS::Traits::GridView,
                                          GFS::Traits::GridView::dimension>
        mapper(gfs.gridView());
    using Dune::PDELab::Backend::native;
    Dune::PDELab::Backend::Vector<GFS, double> v(gfs, 0.0);
    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    LFS lfs(gfs);
    using IC = Dune::PDELab::LFSIndexCache<LFS>;
    IC cache(lfs);
    const auto& ci = st.cutCellInformation();
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> vertexDomainToLinear;
    std::vector<std::set<std::size_t>> cutcellIndexToLinearVertices;
    for (const auto& e : Dune::elements(gfs.gridView())) {
      lfs.bind(e);
      cache.update();
      std::map<std::size_t, std::size_t> domainToCutCellIndex;
      for (unsigned int i = 0; i < cache.size(); ++i) {
        domainToCutCellIndex[cache.dofIndex(i).treeIndex()[1]] = cache.containerIndex(i)[1];
      }
      const auto& geo = e.geometry();
      for (const auto& it : domainToCutCellIndex) {
        for (int i = 0; i < geo.corners(); ++i) {
          auto globalVertexIndex = mapper.subIndex(e, i, GFS::Traits::GridView::dimension);
          auto res = vertexDomainToLinear.insert(
              {{globalVertexIndex, it.first}, vertexDomainToLinear.size()});
          auto cutCellIndex = it.second;
          if (cutCellIndex >= cutcellIndexToLinearVertices.size()) {
            cutcellIndexToLinearVertices.resize(cutCellIndex + 1);
          }
          cutcellIndexToLinearVertices[cutCellIndex].insert(res.first->second);
        }
      }
    }
    using PMatrix =
        Dune::BCRSMatrix<Dune::FieldMatrix<double, Dune::QkStuff::QkSize<1, GFS::Traits::GridView::
                                                                                dimension>::value,
                                           1>>;
    auto matrix = std::make_unique<PMatrix>(cutcellIndexToLinearVertices.size(),
                                                  vertexDomainToLinear.size(), PMatrix::row_wise);
    unsigned int cutCellIndex = 0;
    for (auto it = matrix->createbegin(); it != matrix->createend(); ++it, ++cutCellIndex) {
      for (auto v : cutcellIndexToLinearVertices[cutCellIndex]) {
        it.insert(v);
      }
    }
    for (const auto& e : Dune::elements(gfs.gridView())) {
      lfs.bind(e);
      cache.update();
      auto visitor = cg2udgdetail::make_compute_cg_2_udg_visitor(
          *matrix, e, cache, ci, gfs.gridView(), vertexDomainToLinear);
      Dune::TypeTree::applyToTree(lfs, *visitor);
    }
    return matrix;
  }
}

#endif // DUNEURO_CG_TO_UDG_PROLONGATION_HH
