// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_GEOMETRYADAPTION_HH
#define DUNEURO_GEOMETRYADAPTION_HH

#include <algorithm>
#include <memory>

#include <dune/common/parametertree.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/yaspgrid.hh>
// debug: remove me!
#include <dune/grid/io/file/vtk.hh>

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/volume_conductor.hh>
#if HAVE_NIFTI
#include <duneuro/io/nifti_image_reader.hh>
#endif
#include <duneuro/io/gmsh_tensor_reader.hh>
namespace duneuro
{
  template <class GridView>
  class DeformationFunction
      : public Dune::DiscreteCoordFunction<typename GridView::ctype, GridView::dimension,
                                           DeformationFunction<GridView>>
  {
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using This = DeformationFunction<GridView>;
    using Base = Dune::DiscreteCoordFunction<ctype, dim, This>;

  public:
    DeformationFunction(
        const GridView& gridView,
        std::shared_ptr<std::vector<Dune::FieldVector<ctype, dim>>> deformedPosition)
        : gridView_(gridView), deformedPosition_(deformedPosition)
    {
      assert(deformedPosition_);
      if (deformedPosition_->size() != gridView_.indexSet().size(dim)) {
        DUNE_THROW(Dune::Exception, "every vertex of the grid needs to be deformed");
      }
    }

    void evaluate(const typename GridView::template Codim<dim>::Entity& hostEntity,
                  unsigned int corner, Dune::FieldVector<double, dim>& y) const
    {
      y = (*deformedPosition_)[gridView_.indexSet().index(hostEntity)];
    }

    void evaluate(const typename GridView::template Codim<0>::Entity& hostEntity,
                  unsigned int corner, Dune::FieldVector<double, dim>& y) const
    {
      y = (*deformedPosition_)[gridView_.indexSet().subIndex(hostEntity, corner, dim)];
    }

  private:
    GridView gridView_;

    std::shared_ptr<std::vector<Dune::FieldVector<ctype, dim>>> deformedPosition_;
  };

  template <class GV>
  class VertexToElementsMap
  {
  public:
    using const_iterator = std::vector<std::size_t>::const_iterator;
    enum { dim = GV::dimension };
    using Vertex = typename GV::template Codim<dim>::Entity;

    explicit VertexToElementsMap(const GV& gridView)
        : gridView_(gridView), vertexToElements_(gridView.size(dim))
    {
      Dune::MultipleCodimMultipleGeomTypeMapper<GV>
          elementMapper(gridView, Dune::mcmgElementLayout());
      for (const auto& element : Dune::elements(gridView)) {
        auto elementIndex = elementMapper.index(element);
        for (unsigned int i = 0; i < element.subEntities(dim); ++i) {
          vertexToElements_[gridView.indexSet().subIndex(element, i, dim)].push_back(elementIndex);
        }
      }
    }

    const std::vector<std::size_t>& elements(std::size_t vertexIndex) const
    {
      assert(vertexIndex < vertexToElements_.size());
      return vertexToElements_[vertexIndex];
    }

    const std::vector<std::size_t>& elements(const Vertex& vertex) const
    {
      return elements(gridView_.indexSet().index(vertex));
    }

  private:
    GV gridView_;
    std::vector<std::vector<std::size_t>> vertexToElements_;
  };

  template <class GV>
  class ElementCenters
  {
  public:
    using ctype = typename GV::ctype;
    enum { dim = GV::dimension };
    using Element = typename GV::template Codim<0>::Entity;

    explicit ElementCenters(const GV& gridView)
        : gridView_(gridView), elementToCenter_(gridView_.size(0)),
          elementMapper_(gridView_, Dune::mcmgElementLayout())
    {
      for (const auto& element : Dune::elements(gridView_)) {
        elementToCenter_[elementMapper_.index(element)] = element.geometry().center();
      }
    }

    const Dune::FieldVector<ctype, dim>& center(const Element& element) const
    {
      return center(elementMapper_.index(element));
    }

    const Dune::FieldVector<ctype, dim>& center(std::size_t elementIndex) const
    {
      assert(elementIndex < elementToCenter_.size());
      return elementToCenter_[elementIndex];
    }

    Dune::FieldVector<ctype, dim> centroid(const std::vector<std::size_t>& elements) const
    {
      Dune::FieldVector<ctype, dim> result;
      for (std::size_t e : elements)
        result += center(e);
      result /= elements.size();
      return result;
    }

  private:
    GV gridView_;
    std::vector<Dune::FieldVector<ctype, dim>> elementToCenter_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
  };

  template <class L>
  class LabeledElementSet
  {
  public:
    LabeledElementSet(const std::vector<std::size_t>* elements, const std::vector<L>* labels)
        : elements_(elements), labels_(labels)
    {
      assert(elements);
      assert(labels);
      assert(elements_->size() == labels_->size());
    }

    L minorityLabel() const
    {
      std::map<L, unsigned int> counts;
      for (auto l : *labels_) {
        ++counts[l];
      }
      return std::min_element(
                 counts.begin(), counts.end(),
                 [](const std::pair<L, unsigned int>& l, const std::pair<L, unsigned int>& r) {
                   return l.second < r.second;
                 })
          ->first;
    }

    std::vector<std::size_t> extractElementsWithLabel(L label) const
    {
      std::vector<std::size_t> elements;
      for (unsigned int i = 0; i < elements_->size(); ++i) {
        if ((*labels_)[i] == label) {
          elements.push_back((*elements_)[i]);
        }
      }
      return elements;
    }

    bool hasUniqueMinority() const
    {
      std::map<L, unsigned int> counts;
      for (auto l : *labels_) {
        ++counts[l];
      }
      auto minimal_count =
          std::min_element(counts.begin(), counts.end(), [](const std::pair<L, unsigned int>& l,
                                                            const std::pair<L, unsigned int>& r) {
            return l.second < r.second;
          })->second;
      return std::count_if(
                 counts.begin(), counts.end(),
                 [&](const std::pair<L, unsigned int>& l) { return l.second == minimal_count; })
             == 1;
    }

    std::vector<std::size_t> extractElementsWithMinorityLabel() const
    {
      return extractElementsWithLabel(minorityLabel());
    }

  private:
    const std::vector<std::size_t>* elements_;
    const std::vector<L>* labels_;
  };

  template <class T>
  std::vector<T> extract(const std::vector<T>& values, const std::vector<std::size_t>& indices)
  {
    std::vector<T> result;
    for (const auto& index : indices) {
      if (index >= values.size()) {
        DUNE_THROW(Dune::Exception, "index " << index << " out of bounds (" << values.size()
                                             << ")");
      }
      result.push_back(values[index]);
    }
    return result;
  }

  template <class GV, class L>
  std::shared_ptr<std::vector<Dune::FieldVector<typename GV::ctype, GV::dimension>>>
  compute_deformed_positions(const GV& gridView, const std::vector<L>& elementLabels,
                             const Dune::ParameterTree& config)
  {
    auto deformedPositions =
        std::make_shared<std::vector<Dune::FieldVector<typename GV::ctype, GV::dimension>>>(
            gridView.size(GV::dimension));

    VertexToElementsMap<GV> vertexToElements(gridView);
    ElementCenters<GV> centers(gridView);
    auto shift = config.get<typename GV::ctype>("shift");
    auto minority = config.get<unsigned int>("minority", (1 << (GV::dimension - 1)) - 1);
    for (const auto& vertex : Dune::vertices(gridView)) {
      // extract neighboring elements along with their labels
      const auto& elements = vertexToElements.elements(vertex);
      const auto& labels = extract(elementLabels, elements);
      LabeledElementSet<L> labeledElements(&elements, &labels);
      auto shiftedPosition = vertex.geometry().center();
      if (labeledElements.hasUniqueMinority()) {
        auto minorityElements = labeledElements.extractElementsWithMinorityLabel();
        if (minorityElements.size() <= minority) {
          // compute:  new = old + shift * (centroid - old)
          auto shiftDirection = centers.centroid(minorityElements);
          shiftDirection -= shiftedPosition;
          shiftedPosition.axpy(shift, shiftDirection);
        }
      }
      // store shifted position
      (*deformedPositions)[gridView.indexSet().index(vertex)] = shiftedPosition;
    }
    return deformedPositions;
  }

  template <class HostGrid, class L = std::size_t>
  std::unique_ptr<Dune::GeometryGrid<HostGrid,
                                     DeformationFunction<typename HostGrid::LeafGridView>>>
  create_geometry_adapted_grid(std::shared_ptr<HostGrid> hostGrid,
                               const std::vector<L>& elementLabels,
                               const Dune::ParameterTree& config)
  {
    using HostGridView = typename HostGrid::LeafGridView;
    using CoordFunction = DeformationFunction<HostGridView>;
    using ctype = typename HostGrid::ctype;

    auto hostGridView = hostGrid->leafGridView();
    auto deformedPosition = compute_deformed_positions(hostGridView, elementLabels, config);

    using GridType = Dune::GeometryGrid<HostGrid, CoordFunction>;
    return std::make_unique<GridType>(hostGrid, std::make_shared<CoordFunction>(hostGridView, deformedPosition));
  }

#if HAVE_DUNE_SUBGRID
  template <int dim>
  struct GeometryAdaptedGrid {
    using StructuredGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>;
    using GeometryGrid =
        Dune::GeometryGrid<StructuredGrid,
                           DeformationFunction<typename StructuredGrid::LeafGridView>>;
    using GridType = Dune::SubGrid<dim, GeometryGrid>;

    const GridType& operator*() const
    {
      return *grid;
    }

    const GridType* operator->() const
    {
      return grid.operator->();
    }

    std::unique_ptr<GeometryGrid> geometryGrid;
    std::unique_ptr<GridType> grid;
    std::unique_ptr<std::vector<std::size_t>> labels;
  };

  template <int dim>
  GeometryAdaptedGrid<dim> make_geometry_adapted_grid(const FittedDriverData<dim>& data,
                                                      const Dune::ParameterTree& config)
  {
    auto ygrid = std::make_shared<typename GeometryAdaptedGrid<dim>::StructuredGrid>(
        config.get<Dune::FieldVector<double, dim>>("lower_left"),
        config.get<Dune::FieldVector<double, dim>>("upper_right"),
        config.get<std::array<int, dim>>("cells"));
    auto geometryGrid = create_geometry_adapted_grid(std::move(ygrid), data.labels, config);
    auto ggGv = geometryGrid->leafGridView();
    if (data.labels.size() != static_cast<std::size_t>(ggGv.size(0))) {
      DUNE_THROW(Dune::Exception, "number of labels (" << data.labels.size()
                                                       << ") does not match number of cells ("
                                                       << ggGv.size(0) << ")");
    }
    // subgrid creation
    auto grid = std::make_unique<typename GeometryAdaptedGrid<dim>::GridType>(*geometryGrid);
    grid->createBegin();
    for (const auto& hostElement : Dune::elements(ggGv)) {
      if (data.labels[ggGv.indexSet().index(hostElement)] != 0) {
        grid->insert(hostElement);
      }
    }
    grid->createEnd();
    auto sgGv = grid->leafGridView();
    auto subLabels = std::make_unique<std::vector<std::size_t>>(grid->size(0));
    for (const auto& subElement : Dune::elements(sgGv)) {
      (*subLabels)[sgGv.indexSet().index(subElement)] =
          data.labels[ggGv.indexSet().index(grid->template getHostEntity<0>(subElement))] - 1;
    }
    return GeometryAdaptedGrid<dim>{std::move(geometryGrid), std::move(grid), std::move(subLabels)};
  }

  template <int dim>
  std::shared_ptr<VolumeConductor<typename GeometryAdaptedGrid<dim>::GridType>>
  make_geometry_adapted_volume_conductor(
      std::unique_ptr<typename GeometryAdaptedGrid<dim>::GridType> grid,
      std::unique_ptr<std::vector<std::size_t>> labels, const std::vector<double>& conductivities,
      const Dune::ParameterTree& config)
  {
    using VC = VolumeConductor<typename GeometryAdaptedGrid<dim>::GridType>;
    using Tensor = typename VC::TensorType;
    std::vector<Tensor> tensors;
    for (auto value : conductivities) {
      Tensor t;
      for (unsigned int r = 0; r < t.N(); ++r) {
        for (unsigned int c = 0; c < t.M(); ++c) {
          t[r][c] = r == c ? value : 0.0;
        }
      }
      tensors.push_back(t);
    }
    return std::make_shared<VC>(std::move(grid), *labels, tensors);
  }

#if HAVE_NIFTI
  template <int dim>
  struct GeometryAdaptedGridReader {
    static GeometryAdaptedGrid<dim> read(const Dune::ParameterTree& config)
    {
      std::vector<std::size_t> labels;
      std::array<unsigned int, dim> cells;
      NiftiImageReader::read<dim>(config.get<std::string>("filename"), std::back_inserter(labels),
                                  cells);
      std::array<int, dim> s;
      std::copy(cells.begin(), cells.end(), s.begin());
      auto ygrid = std::make_shared<typename GeometryAdaptedGrid<dim>::StructuredGrid>(
          config.get<Dune::FieldVector<double, dim>>("lower_left"),
          config.get<Dune::FieldVector<double, dim>>("upper_right"), s);
      // geometry adaption
      auto geometryGrid = create_geometry_adapted_grid(std::move(ygrid), labels, config);
      auto ggGv = geometryGrid->leafGridView();
      assert(labels.size() == static_cast<std::size_t>(ggGv.size(0)));
      // subgrid creation
      auto grid =
          std::make_unique<typename GeometryAdaptedGrid<dim>::GridType>(*geometryGrid);
      grid->createBegin();
      for (const auto& hostElement : Dune::elements(ggGv)) {
        if (labels[ggGv.indexSet().index(hostElement)] != 0) {
          grid->insert(hostElement);
        }
      }
      grid->createEnd();
      auto sgGv = grid->leafGridView();
      auto subLabels = std::make_unique<std::vector<std::size_t>>(grid->size(0));
      for (const auto& subElement : Dune::elements(sgGv)) {
        (*subLabels)[sgGv.indexSet().index(subElement)] =
            labels[ggGv.indexSet().index(grid->template getHostEntity<0>(subElement))] - 1;
      }
      return GeometryAdaptedGrid<dim>{std::move(geometryGrid), std::move(grid),
                                      std::move(subLabels)};
    }
  };

#endif
#endif
}

#endif // DUNEURO_GEOMETRYADAPTION_HH
