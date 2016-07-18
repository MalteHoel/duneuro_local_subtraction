#ifndef DUNEURO_VOLUMECONDUCTOR_HH
#define DUNEURO_VOLUMECONDUCTOR_HH

#include <functional>
#include <memory>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/utility/hierarchicsearch.hh>

namespace duneuro
{
  template <class GV, class T, class I = std::size_t>
  class IndirectEntityMapping
  {
  public:
    typedef T TensorType;
    typedef typename GV::template Codim<0>::Entity EntityType;

    IndirectEntityMapping(const GV& gridView, const std::vector<T>& tensors,
                          const std::vector<I>& indexToTensor)
        : gridView_(gridView), mapper_(gridView), tensors_(tensors), indexToTensor_(indexToTensor)
    {
      assert(tensors.size() > 0);
      assert(indexToTensor_.size() == static_cast<std::size_t>(gridView_.size(0)));
    }

    const T& operator()(const EntityType& e) const
    {
      auto index = mapper_.index(e);
      assert(index < indexToTensor_.size());
      auto tensorIndex = indexToTensor_[index];
      assert(tensorIndex < tensors_.size());
      return tensors_[tensorIndex];
    }

  private:
    // copy the grid view so that the mapper is always valid
    GV gridView_;
    Dune::SingleCodimSingleGeomTypeMapper<GV, 0> mapper_;
    std::vector<T> tensors_;
    // maps an entity index to the tensor index
    std::vector<I> indexToTensor_;
  };

  template <class GV, class T>
  class DirectEntityMapping
  {
  public:
    typedef T TensorsType;
    typedef typename GV::template Codim<0>::Entity EntityType;

    DirectEntityMapping(const GV& gridView, const std::vector<T>& tensors)
        : gridView_(gridView), mapper_(gridView_), tensors_(tensors)
    {
      assert(gridView_.size(0) == tensors_.size());
    }

    const T& operator()(const EntityType& e) const
    {
      auto index = mapper_.index(e);
      assert(index < tensors_.size());
      return tensors_[index];
    }

  private:
    // copy the grid view so that the mapper is always valid
    GV gridView_;
    Dune::SingleCodimSingleGeomTypeMapper<GV, 0> mapper_;
    std::vector<T> tensors_;
  };

  template <class G>
  class VolumeConductor
  {
  public:
    typedef G GridType;
    enum { dim = G::dimension };
    typedef typename G::ctype ctype;
    typedef typename G::template Codim<0>::Entity EntityType;
    typedef Dune::FieldMatrix<ctype, dim, dim> TensorType;
    typedef std::function<const TensorType&(const EntityType&)> MappingType;
    typedef typename G::LeafGridView GridView;

    VolumeConductor(std::unique_ptr<G> grid, std::unique_ptr<MappingType> mapper)
        : grid_(std::move(grid)), mapper_(std::move(mapper)), gridView_(grid_->leafGridView())
    {
    }

    const G& grid() const
    {
      return *grid_;
    }

    G& grid()
    {
      return *grid_;
    }

    const GridView& gridView() const
    {
      return gridView_;
    }

    const TensorType& tensor(const EntityType& entity) const
    {
      return (*mapper_)(entity);
    }

    G* releaseGrid()
    {
      return grid_.release();
    }

    MappingType* releaseMapping()
    {
      return mapper_.release();
    }

  private:
    std::unique_ptr<G> grid_;
    std::unique_ptr<MappingType> mapper_;
    GridView gridView_;
  };
}

#endif // DUNEURO_VOLUMECONDUCTOR_HH
