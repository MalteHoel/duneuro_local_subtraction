#ifndef DUNEURO_ELEMENT_PATCH_HH
#define DUNEURO_ELEMENT_PATCH_HH

#include <memory>
#include <vector>

#include <dune/common/parametertree.hh>

#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <duneuro/common/element_neighborhood_map.hh>

namespace duneuro
{
  enum class ElementPatchInitialization { singleElement, closestVertex };
  enum class ElementPatchExtension { vertex, intersection };

  static inline ElementPatchInitialization
  elementPatchInitializationFromString(const std::string& name)
  {
    if (name == "single_element") {
      return ElementPatchInitialization::singleElement;
    } else if (name == "closest_vertex") {
      return ElementPatchInitialization::closestVertex;
    } else {
      DUNE_THROW(Dune::Exception, "unknown element patch initialization \"" << name << "\"");
    }
  }

  static inline ElementPatchExtension elementPatchExtensionFromString(const std::string& name)
  {
    if (name == "vertex") {
      return ElementPatchExtension::vertex;
    } else if (name == "intersection") {
      return ElementPatchExtension::intersection;
    } else {
      DUNE_THROW(Dune::Exception, "unknown element patch extension \"" << name << "\"");
    }
  }

  template <typename GV>
  class ElementPatch
  {
  public:
    using GridView = GV;
    using Coordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using Element = typename GV::template Codim<0>::Entity;
    using Intersection = typename GV::Intersection;
    using ElementMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, 0>;
    using VertexMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, GV::dimension>;

    template <typename ElementSearch>
    ElementPatch(std::shared_ptr<ElementNeighborhoodMap<GV>> elementNeighborhoodMap,
                 const ElementSearch& elementSearch, const Coordinate& position,
                 ElementPatchInitialization initialization,
                 std::function<bool(Element)> elementFilter = [](const Element&) { return true; })
        : elementNeighborhoodMap_(elementNeighborhoodMap)
        , elementFilter_(elementFilter)
        , elementMapper_(elementNeighborhoodMap_->gridView())
        , vertexMapper_(elementNeighborhoodMap_->gridView())
    {
      switch (initialization) {
      case ElementPatchInitialization::singleElement:
        initializeSingleElement(elementSearch, position);
        break;
      case ElementPatchInitialization::closestVertex:
        initializeClosestVertex(elementSearch, position);
        break;
      }
    }

    void extend(ElementPatchExtension extension)
    {
      std::vector<Element> candidates;
      switch (extension) {
      case ElementPatchExtension::vertex:
        for (const auto& element : elements_) {
          elementNeighborhoodMap_->getVertexNeighbors(element, std::back_inserter(candidates));
        }
        break;
      case ElementPatchExtension::intersection:
        for (const auto& element : elements_) {
          elementNeighborhoodMap_->getIntersectionNeighbors(element,
                                                            std::back_inserter(candidates));
        }
        break;
      }
      for (const auto& candidate : candidates) {
        auto index = elementMapper_.index(candidate);
        if ((elementIndices_.count(index) == 0) && elementFilter_(candidate)) {
          elements_.push_back(candidate);
          elementIndices_.insert(index);
        }
      }
    }

    void extend(const std::vector<ElementPatchExtension>& extensions, std::size_t repeatUntil = 0)
    {
      auto old = elements_.size();
      for (const auto& type : extensions) {
        extend(type);
      }
      while (elements_.size() < repeatUntil && old != elements_.size()) {
        old = elements_.size();
        for (const auto& type : extensions) {
          extend(type);
        }
      }
    }

    void extend(const std::vector<std::string>& extensions, std::size_t repeatUntil = 0)
    {
      std::vector<ElementPatchExtension> ex;
      std::transform(extensions.begin(), extensions.end(), std::back_inserter(ex),
                     elementPatchExtensionFromString);
      extend(ex, repeatUntil);
    }

    //
    // I deleted the reference here
    //
    const std::vector<Element> elements() const
    {
      return elements_;
    }

    bool contains(const Element& element) const
    {
      return elementIndices_.count(elementMapper_.index(element)) > 0;
    }

    template <class I>
    void extractBoundaryIntersections(I out) const
    {
      for (const auto& element : elements_) {
        for (const auto& is : Dune::intersections(elementNeighborhoodMap_->gridView(), element)) {
          if (is.neighbor() && !contains(is.outside())) {
            *out++ = is;
          }
        }
      }
    }

    std::vector<Intersection> extractBoundaryIntersections() const
    {
      std::vector<Intersection> out;
      extractBoundaryIntersections(std::back_inserter(out));
      return out;
    }


    // written by me
    Element initialElement() 
    {
      return element_of_start_position_;
    }
    
    size_t sizeOfPatch() {
      return elements_.size();
    }

    // For the CG localized subtraction source model, we need one additional vertex extension as a transitional region 
    // for chi to drop to zero. This function computes a vector containing the elements of this transitional region,
    // where we assume that the inner region of the patch has already been computed
    std::vector<Element> computeTransitionRegion() 
    {
      // we essentially need to perfrom one additional vertex extension and simply store the new elements inside the array
      
      // we first get all elements that share a vertex with one element in the current patch. Note that candidates get included mutliple times
      std::vector<Element> candidates;
      for(const auto& element : elements_) {
      	elementNeighborhoodMap_->getVertexNeighbors(element, std::back_inserter(candidates));
      }
      
      // we now check all candidates to see if they are new additions
      std::vector<Element> transitionElements;
      std::set<std::size_t> visitedTransitionElementIndices;
      for(const auto& candidate : candidates) {
        auto index = elementMapper_.index(candidate);
        // check if candidate is not in inner region as was not already included earlier
        if(elementIndices_.count(index) == 0 && visitedTransitionElementIndices.count(index) == 0) {
          transitionElements.push_back(candidate);
          visitedTransitionElementIndices.insert(index);
        }
      }
      
      return transitionElements;  
    }

    // end of part written by me
  
  private:
    std::shared_ptr<ElementNeighborhoodMap<GV>> elementNeighborhoodMap_;
    std::function<bool(Element)> elementFilter_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    // written by me
    Element element_of_start_position_;

    std::vector<Element> elements_;
    std::set<std::size_t> elementIndices_;

    template <typename ElementSearch>
    void initializeSingleElement(const ElementSearch& elementSearch, const Coordinate& position)
    {
      element_of_start_position_ = elementSearch.findEntity(position);
      if (elementFilter_(element_of_start_position_)) {
        elements_.push_back(element_of_start_position_);
        elementIndices_.insert(elementMapper_.index(element_of_start_position_));
      }
    }

    template <typename ElementSearch>
    void initializeClosestVertex(const ElementSearch& elementSearch, const Coordinate& position)
    {
      element_of_start_position_ = elementSearch.findEntity(position);
      const auto& geo = element_of_start_position_.geometry();
      // find closest corner
      unsigned int minCorner = 0;
      double minDistance = std::numeric_limits<double>::max();
      for (unsigned int i = 0; i < geo.corners(); ++i) {
        Coordinate tmp = position;
        tmp -= geo.corner(i);
        double tn = tmp.two_norm();
        if (tn < minDistance) {
          minDistance = tn;
          minCorner = i;
        }
      }
      // retrieve elements belonging to that corner
      std::vector<Element> candidates;
      elementNeighborhoodMap_->getNeighborsOfVertex(
          vertexMapper_.subIndex(element_of_start_position_, minCorner, GV::dimension),
          std::back_inserter(candidates));
      // filter them and push them to the list
      for (const auto& e : candidates) {
        if (elementFilter_(e)) {
          elements_.push_back(e);
          elementIndices_.insert(elementMapper_.index(e));
        }
      }
    }
  };

  template <class VC, class ES>
  std::function<bool(typename VC::EntityType)>
  make_element_filter(std::shared_ptr<const VC> volumeConductor, const ES& elementSearch,
                      const Dune::FieldVector<typename VC::ctype, VC::dim>& position, bool restrict)
  {
    if (restrict) {
      auto reference = volumeConductor->tensor(elementSearch.findEntity(position));
      return [volumeConductor, reference](const typename VC::EntityType& e) {
        auto diff = volumeConductor->tensor(e);
        diff -= reference;
        return diff.frobenius_norm2() < 1e-8;
      };
    } else {
      return [](const typename VC::EntityType&) { return true; };
    }
  }

  template <class VC, class ES>
  std::unique_ptr<ElementPatch<typename VC::GridView>> make_element_patch(
      std::shared_ptr<const VC> volumeConductor,
      std::shared_ptr<ElementNeighborhoodMap<typename VC::GridView>> elementNeighborhoodMap,
      const ES& elementSearch, const Dune::FieldVector<typename VC::ctype, VC::dim>& position,
      const Dune::ParameterTree& config)
  {
    auto patch = std::make_unique<ElementPatch<typename VC::GridView>>(
        elementNeighborhoodMap, elementSearch, position,
        elementPatchInitializationFromString(config.get<std::string>("initialization")),
        make_element_filter(volumeConductor, elementSearch, position,
                            config.get<bool>("restrict")));
    patch->extend(config.get("extensions", std::vector<std::string>()),
                  config.get<unsigned int>("repeat_until", 0));
    return patch;
  }

#if HAVE_DUNE_UDG
  template <class ST, class ES>
  std::function<bool(typename ST::Entity)>
  make_element_filter(std::shared_ptr<ST> subTriangulation, const ES& elementSearch,
                      const Dune::FieldVector<typename ST::ctype, ST::dim>& position,
                      unsigned int domainIndex)
  {
    return [subTriangulation, domainIndex](const typename ST::Entity& e) {
      return subTriangulation->isHostCell(e, domainIndex);
    };
  }

  template <class ST, class ES>
  std::unique_ptr<ElementPatch<typename ST::GridView>> make_element_patch(
      std::shared_ptr<ST> subTriangulation,
      std::shared_ptr<ElementNeighborhoodMap<typename ST::GridView>> elementNeighborhoodMap,
      const ES& elementSearch,
      const Dune::FieldVector<typename ST::GridView::ctype, ST::GridView::dimension>& position,
      unsigned int positionDomain, const Dune::ParameterTree& config)
  {
    auto patch = std::make_unique<ElementPatch<typename ST::GridView>>(
        elementNeighborhoodMap, elementSearch, position,
        elementPatchInitializationFromString(config.get<std::string>("initialization")),
        make_element_filter(subTriangulation, elementSearch, position, positionDomain));
    patch->extend(config.get("extensions", std::vector<std::string>()),
                  config.get<unsigned int>("repeat_until", 0));
    return patch;
  }
#endif

  template <class EP>
  class ElementPatchVTKFunction : public Dune::VTKFunction<typename EP::GridView>
  {
    using GV = typename EP::GridView;
    typedef typename GV::ctype DF;
    enum { n = GV::dimension };
    typedef typename GV::template Codim<0>::Entity Entity;

  public:
    ElementPatchVTKFunction(const EP& ep, std::string name) : ep_(ep), name_(name)
    {
    }

    virtual int ncomps() const
    {
      return 1;
    }

    virtual double evaluate(int comp, const Entity& e, const Dune::FieldVector<DF, n>& xi) const
    {
      return ep_.contains(e) ? 1 : 0;
    }

    virtual std::string name() const
    {
      return name_;
    }

  private:
    const EP& ep_;
    std::string name_;
  };

  template <class EP>
  std::unique_ptr<ElementPatchVTKFunction<EP>>
  make_element_patch_vtk_function(const EP& patch, const std::string& name = "patch")
  {
    return std::make_unique<ElementPatchVTKFunction<EP>>(patch, name);
  }
}

#endif // DUNEURO_ELEMENT_PATCH_HH
