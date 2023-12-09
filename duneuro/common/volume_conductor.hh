#ifndef DUNEURO_VOLUMECONDUCTOR_HH
#define DUNEURO_VOLUMECONDUCTOR_HH

#include <functional>
#include <memory>
#include <vector>
#include <set>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <duneuro/common/element_neighborhood_map.hh>

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
    typedef typename G::template Codim<dim>::Entity VertexType;
    typedef Dune::FieldMatrix<ctype, dim, dim> TensorType;
    typedef typename G::LeafGridView GridView;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView, dim>::Index VertexIndex;

    VolumeConductor(std::unique_ptr<G> grid, std::vector<std::size_t> labels,
                    std::vector<TensorType> tensors)
        : grid_(std::move(grid))
        , labels_(labels)
        , tensors_(tensors)
        , gridView_(grid_->leafGridView())
        , elementMapper_(gridView_)
        , vertexMapper_(gridView_)
        , elementNeighborhoodMapPtr_(nullptr)
        , elementNeighborhoodMapComputed_(false)
        , hasInsertionIndices_(false)
        , elementInsertionIndices_(0)
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
    
    // if VolumeConductor was constructed by directly specifying nodes and elements, we want to keep track of the insertion indices
    VolumeConductor(std::unique_ptr<G> grid, std::vector<std::size_t> labels,
                    std::vector<TensorType> tensors, std::vector<std::size_t> elementInsertionIndices, std::vector<std::size_t> vertexInsertionIndices)
        : grid_(std::move(grid))
        , labels_(labels)
        , tensors_(tensors)
        , gridView_(grid_->leafGridView())
        , elementMapper_(gridView_)
        , vertexMapper_(gridView_)
        , elementNeighborhoodMapPtr_(nullptr)
        , elementNeighborhoodMapComputed_(false)
        , hasInsertionIndices_(true)
        , elementInsertionIndices_(elementInsertionIndices)
        , vertexInsertionIndices_(vertexInsertionIndices)
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
    
    std::size_t insertionIndex(const EntityType& entity) const
    {
      if(hasInsertionIndices_) {
        return elementInsertionIndices_[elementMapper_.index(entity)];
      }
      else {
        DUNE_THROW(Dune::Exception, "this volume conductor does not store insertion indices");
      } 
    }
    
    std::size_t vertexInsertionIndex(const VertexType& vertex) const
    {
      if(hasInsertionIndices_) {
        return vertexInsertionIndices_[vertexMapper_.index(vertex)];
      }
      else {
        DUNE_THROW(Dune::Exception, "this volume conductor does not store insertion indices");
      } 
    }

    G* releaseGrid()
    {
      return grid_.release();
    }

    void computeElementNeighborhoodMap()
    {
      elementNeighborhoodMapPtr_ = std::make_shared<ElementNeighborhoodMap<GridView>>(gridView_);
      elementNeighborhoodMapComputed_ = true;
    }

    std::shared_ptr<ElementNeighborhoodMap<GridView>> elementNeighborhoodMap() const
    {
      if(elementNeighborhoodMapComputed_) {
        return elementNeighborhoodMapPtr_;
      }
      else {
        DUNE_THROW(Dune::Exception, "Element neighborhood map needed, but not computed");
      }
    }
    
    std::set<VertexIndex> venantVertices(const std::set<std::size_t>& sourceCompartments) const 
    {
      const auto& vertexToElements = elementNeighborhoodMapPtr_->vertexToElements();
      std::set<VertexIndex> venantVertexSet;
      
      for(const auto& vertex : vertices(gridView_)) {
        bool venantVertex = true;
        for(auto elementSeed : vertexToElements[vertexMapper_.index(vertex)]) {
          if(sourceCompartments.find(labels_[elementMapper_.index(grid_->entity(elementSeed))]) == sourceCompartments.end()) {
            venantVertex = false;
            break;
          }
        }
        if(venantVertex) {
          venantVertexSet.insert(vertexMapper_.index(vertex));
        }
      }
      
      return venantVertexSet;
    }
    
    const std::vector<TensorType>& tensors() const {
      return tensors_;
    }
    
    VertexIndex vertexIndex(typename GridView::template Codim<dim>::Entity::EntitySeed vertexSeed) const
    {
      return vertexMapper_.index(grid_->entity(vertexSeed));
    }

  private:
    std::unique_ptr<G> grid_;
    const std::vector<std::size_t> labels_;
    const std::vector<TensorType> tensors_;
    GridView gridView_;
    Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> elementMapper_;
    Dune::SingleCodimSingleGeomTypeMapper<GridView, dim> vertexMapper_;
    std::shared_ptr<ElementNeighborhoodMap<GridView>> elementNeighborhoodMapPtr_;
    bool elementNeighborhoodMapComputed_;
    bool hasInsertionIndices_;
    std::vector<std::size_t> elementInsertionIndices_;
    std::vector<std::size_t> vertexInsertionIndices_;
  };
}

#endif // DUNEURO_VOLUMECONDUCTOR_HH
