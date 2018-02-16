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
  template <class G>
  class VolumeConductor
  {
  public:
    typedef G GridType;
    enum { dim = G::dimension };
    typedef typename G::ctype ctype;
    typedef typename G::template Codim<0>::Entity EntityType;
    typedef Dune::FieldMatrix<ctype, dim, dim> TensorType;
    typedef typename G::LeafGridView GridView;

    VolumeConductor(std::unique_ptr<G> grid, std::vector<std::size_t> labels,
                    std::vector<TensorType> tensors)
        : grid_(std::move(grid))
        , labels_(labels)
        , tensors_(tensors)
        , gridView_(grid_->leafGridView())
        , elementMapper_(gridView_)
    {
      // check if we are given one label for each element
      if (labels.size() != elementMapper_.size()) {
        DUNE_THROW(Dune::Exception, "number of labels ("
                                        << labels.size()
                                        << ") has to match the number of elements ("
                                        << elementMapper_.size() << ")");
      }
      // check if all labels are valid
      for (const auto& l : labels) {
        if (l >= tensors_.size()) {
          DUNE_THROW(Dune::Exception,
                     "label " << l << " is not a valid index into the tensors vector of size "
                              << tensors_.size());
        }
      }
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
      return tensors_[label(entity)];
    }

    std::size_t label(const EntityType& entity) const
    {
      return labels_[elementMapper_.index(entity)];
    }

    G* releaseGrid()
    {
      return grid_.release();
    }

  private:
    std::unique_ptr<G> grid_;
    const std::vector<std::size_t> labels_;
    const std::vector<TensorType> tensors_;
    GridView gridView_;
    Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> elementMapper_;
  };
}

#endif // DUNEURO_VOLUMECONDUCTOR_HH
