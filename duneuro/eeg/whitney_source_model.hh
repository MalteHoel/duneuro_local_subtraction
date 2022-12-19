#ifndef DUNEURO_WHITNEY_SOURCE_MODEL_HH
#define DUNEURO_WHITNEY_SOURCE_MODEL_HH

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <array>
#include <vector>
#include <queue>
#include <functional>
#include <utility>
#include <set>
#include <cmath>
#include <limits>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/eeg/source_model_interface.hh>

///////////////////////////////////////////////////////////////////////
//                                                                   //
// WhitneySourceModel creates a source model with "Whitney-          //
// elements, taking a combination of face intersecting and edgewise  //
// dipolar sources and estimating the original dipole position and   //
// moment by creating a linear combination of these source dipoles   //
// by interpolating. The script follows the methods presented in     //
// Pursiainen, S., Vorwek, J. and Wolters, C.,                       //
// "Electroencephalography (EEG) forward modeling                    //
// via H(div) Finite Element Sources With Focal Interpolation",      //
// published in Physic in Medicine and Biology, vol 61, no. 24,      //
// pp. 8502-8520, 11 2016                                            //
// and                                                               //
// Miinalainen, T., Rezaei, A., Us, D., Nuessing, A., Engwer,        //
// C., Wolters, C., Pursiainen, S. "A realistic, accurate and fast   //
// source modeling approach for the EEG forward problem", published  //
// in Neuroimage, vol 184, pp. 56-67, 01 2019                        //
//                                                                   //
///////////////////////////////////////////////////////////////////////


// Required parameters:
// std::shared_ptr<VC> volumeConductor   - Volume Conductor
// const GFS gfs  - Global function space
// mutable CacheType cache   - cache element
// const Real referenceLength - Reference length for PBO interpolation
// (usually double  - three times) the length of the longest edge in the mesh
// const bool restricted - boolean if restricting the sources in one layer
// const WhitneyFaceBased - variable defining if element's face based dipoles
// are included:
// "all" is 4 dipoles, or "none" 0
// const WhitneyEdgeBased - variable defining if element's or neighboring edge
// based dipoles are included:
// "all" takes neighboring and elements edge dipoles, "internal" just element's
// edges and "none" is for 0
// const std::string interpolation - name of the interpolation, either
// 'MPO' or 'PBO'

// Currently source model functions only with tetrahedras. Also, it is
// required that the basis functions are nodal basis functions

namespace duneuro {

// the central idea in the Whitney approach consists in approximating a given dipole
// M \delta_{x_0} in the form M \delta_{x_0} \approx \sum_{l = 1}^k c_l w_l, where
// w_1, ..., w_k Elements of H(div). On then uses the sum on the right hand side for
// the finite element discretization. In the interpolation and assembly methods we want
// to use duck typing for the functions w_1, ..., w_k. The following describes the
// interface we expect an implementation to adhere to.
/*
class HdivFunction {
public:

  // Static information exported by this class
  using DipoleType = ...;
  using Scalar = ...;
  enum {dim = ...};

  // the interpolation methods described in the papers cited above rely on associating
  // a synthetic dipole to each given H(div) basis functions. More concretely, they interpolate
  // a given dipole by finding a "similar" linear combination of the synthetic dipoles, where
  // similar depends on the interpolation approach
  DipoleType syntheticDipole() const;

  // Once we have an approximation M \delta_{x_0} = \sum_{l = 1}^k c_l w_l, we have to compute
  // \int_{\Omega} \sum_{l = 1} c_l w_l \psi_i dV for every test function \psi_i. If w is the
  // function described by an instance of this class, this method is supposed to perform
  // vector[i] += \int_{\Omega} coefficient * w * \psi_i
  // for every test function \psi_i.
  template<class DOFVector>
  void accumulate(DOFVector& output, const Scalar coefficient) const;
}; // class HdivFunction
*/


// In the papers cited above so called "face intersecting" (FI) and "edgewise" (EW) divergence
// conforming functions are used. It turns out that all necessary informations for both the FI
// and the EW functions can be deduced from specifying two vertices in the grid. For the FI
// functions these are the two nodes opposite to the given triangular face, while for the EW
// functions these are the nodes making up the edge
template<class GridFunctionSpace>
class TwoVerticesHdivFunction
{
public:
  using GridView = typename GridFunctionSpace::Traits::GridViewType;
  enum {dim = GridView::dimension};
  using Scalar = typename GridView::ctype;
  using Coordinate = Dune::FieldVector<Scalar, dim>;
  using DipoleType = Dipole<Scalar, dim>;
  using VertexEntity = typename GridView::template Codim<dim>::Entity;
  using VertexSeed = typename GridView::Grid::template Codim<dim>::EntitySeed;
  using Cache = typename Dune::PDELab::EntityIndexCache<GridFunctionSpace>;

  TwoVerticesHdivFunction(const GridFunctionSpace& gfs, const VertexEntity& vertex_entity_0, const VertexEntity& vertex_entity_1)
    : gfs_(gfs), vertex_seeds_({vertex_entity_0.seed(), vertex_entity_1.seed()})
    , vertex_0_(vertex_entity_0.geometry().center())
    , vertex_1_(vertex_entity_1.geometry().center())
    , normed_diff_(vertex_1_ - vertex_0_)
    , distance_(normed_diff_.two_norm())
  {
    normed_diff_ /= distance_;
  }

  DipoleType syntheticDipole() const {
    return Dipole(0.5 * (vertex_0_ + vertex_1_), normed_diff_);
  }

  template <class DOFVector>
  void accumulate(DOFVector& vec, const Scalar coefficient) const
  {
    Cache cache(gfs_);

    const auto& grid = gfs_.gridView().grid();
    auto vertex_entity_0 = grid.entity(vertex_seeds_[0]);
    auto vertex_entity_1 = grid.entity(vertex_seeds_[1]);

    Scalar update = coefficient / distance_;

    cache.update(vertex_entity_0);
    vec[cache.containerIndex(0)] -= update;

    cache.update(vertex_entity_1);
    vec[cache.containerIndex(0)] += update;
  }

private:
  std::array<VertexSeed, 2> vertex_seeds_;
  const GridFunctionSpace& gfs_;
  Coordinate vertex_0_;
  Coordinate vertex_1_;
  Coordinate normed_diff_;
  Scalar distance_;
}; // class TwoVerticesHdivFunction

template<class BasisFunction>
class SourceConfigurationBuilderInterface
{
public:

  virtual std::vector<BasisFunction> build(const typename BasisFunction::GridView::template Codim<0>::Entity& initialElement,
                                           const typename BasisFunction::DipoleType& dipole) const = 0;

  virtual ~SourceConfigurationBuilderInterface() {}
}; // SourceConfigurationBuilderInterface

template<class VolumeConductor, class GridFunctionSpace>
class FixedSourceConfigurationBuilder
  : public SourceConfigurationBuilderInterface<TwoVerticesHdivFunction<GridFunctionSpace>>
{
public:
  using BasisFunction =  TwoVerticesHdivFunction<GridFunctionSpace>;
  using DipoleType = typename BasisFunction::DipoleType;
  using ElementType = typename GridFunctionSpace::Traits::GridView::template Codim<0>::Entity;
  enum {dim = GridFunctionSpace::Traits::GridViewType::dimension};
  enum class WhitneyEdges {none, internal, all};
  enum {NUMBER_OF_EDGES = 6};
  enum {NUMBER_OF_EDGES_WITH_NEIGHBORS = 18};
  enum {NUMBER_OF_FACES = 4};

  FixedSourceConfigurationBuilder(std::shared_ptr<const VolumeConductor> volumeConductor,
                                  const GridFunctionSpace& gfs,
                                  const Dune::ParameterTree& config)
    : volumeConductor_(volumeConductor)
    , gfs_(gfs)
    , restricted_(config.get<bool>("restricted"))
    , whitneyFaces_(config.get<bool>("faceSources"))
    , whitneyEdges_(WhitneyEdgesFromString(config.get<std::string>("edgeSources")))
    , basisSizeEstimate_(0)
  {
    if (whitneyFaces_) {
      basisSizeEstimate_ += NUMBER_OF_FACES;
    }
    if (whitneyEdges_ == WhitneyEdges::internal) {
      basisSizeEstimate_ += NUMBER_OF_EDGES;
    }
    else if (whitneyEdges_ == WhitneyEdges::all) {
      basisSizeEstimate_ += NUMBER_OF_EDGES_WITH_NEIGHBORS;
    }
  }

  // The dipole argument is not used, but is needed to adhere to the interface
  std::vector<BasisFunction> build(const ElementType& initialElement, const DipoleType& dipole) const override
  {
    std::vector<BasisFunction> basisFunctions;
    basisFunctions.reserve(basisSizeEstimate_);

    auto ref = referenceElement(initialElement.geometry());

    // first add internal edges, if required
    if(whitneyEdges_ == WhitneyEdges::internal || whitneyEdges_ == WhitneyEdges::all) {
      for(int edge_index = 0; edge_index < ref.size(dim - 1); ++edge_index) {
        auto vertexIndexRange = ref.subEntities(edge_index, dim - 1, dim);
        auto vertexIndexIterator = vertexIndexRange.begin();
        auto index_1 = *(vertexIndexIterator++);
        auto index_2 = *(vertexIndexIterator);
        basisFunctions.emplace_back(gfs_, initialElement.template subEntity<dim>(index_1), initialElement.template subEntity<dim>(index_2));
      }
    }

    // now potentially add basis functions from adjacent elements
    if(whitneyFaces_ || whitneyEdges_ == WhitneyEdges::all) {
      auto initialConductivity = volumeConductor_->tensor(initialElement);
      for(const auto& intersection : intersections(gfs_.gridView(), initialElement)) {
        if(!intersection.neighbor()) continue;

        auto outsideEntity = intersection.outside();
        if(restricted_ && volumeConductor_->tensor(outsideEntity) != initialConductivity) continue;

        // the following assumes tetrahedral elements. It just so happens that the element numberings
        // of the Dune tetrahedra have the property that if i is the local index of a vertex
        // and j is the local index of the face opposite to this vertex, we have i + j = 3.
        auto outsideVertexIndex = dim - intersection.indexInOutside();

        if(whitneyFaces_) {
          auto insideVertexIndex = dim - intersection.indexInInside();
          basisFunctions.emplace_back(gfs_, initialElement.template subEntity<dim>(insideVertexIndex), outsideEntity.template subEntity<dim>(outsideVertexIndex));
        }

        if(whitneyEdges_ == WhitneyEdges::all) {
          for(const auto& vertexIndex : ref.subEntities(intersection.indexInOutside(), 1, dim)) {
            basisFunctions.emplace_back(gfs_, outsideEntity.template subEntity<dim>(vertexIndex), outsideEntity.template subEntity<dim>(outsideVertexIndex));
          }
        }
      }
    } // potentialled added basis functions from adjacent elements

    return basisFunctions;
  }

private:
  WhitneyEdges WhitneyEdgesFromString(const std::string& value)
  {
    if(value == "none")
      return WhitneyEdges::none;
    else if (value == "internal")
      return WhitneyEdges::internal;
    else if (value == "all")
      return WhitneyEdges::all;
    else
      DUNE_THROW(Dune::Exception, "edgeSources has invalid value " << value << ", please choose one from {none, internal, all}");
  }

  std::shared_ptr<const VolumeConductor> volumeConductor_;
  const GridFunctionSpace& gfs_;
  bool restricted_;
  bool whitneyFaces_;
  WhitneyEdges whitneyEdges_;
  size_t basisSizeEstimate_;
}; // FixedSourceConfigurationBuilder

// The new idea in Miinalainen et. al., 2018 is to dynamically construct the source configuration.
// The idea is as follows: One supplies an additional parameter, called n_elems. We now start from the element
// containing the dipole, add all internal edge sources. This is the first element. We then iteratively perform the
// following steps.
// Assume we have a patch consisting of k elements. We then use some heuristic to choose an element which shares at
// least one face with one of the k face elements. We then add this element to the patch, and include all Whitney basis
// functions corresponding to
//    1) an edge of the new element or
//    2) to a face of the new element with one of the k given patch elements
// into the source configuration, given they are not already contained. We repeat this until we have reached n_elems patch elements.
// Now we need to describe an heuristic for choosing the next patch element.
// In the paper, Miinalainen et. al. are quite vague about this. Looking at the description on top of page 4 of that paper,
// and at figures 1 and 2 in the paper, one gets the impression that all facial neighbors of the dipole elmement are choosen,
// as long as they are in the same compartment as the source. As far as I can tell, it is not described how one should proceed if
// one cannot reach the desired amount of elements this way.
// In Tuulis implementation on the other hand, another approach is taken. After the initial element, only elements such that the element itself
// and all of its facial neighbors are inside the same compartment as the dipole are added. This is first tested for all facial neighbors of
// the initial element, and these are added if this condition is fulfilled. If the desired number of elements is not reached after this step,
// Tuulis implementation then loops over all facial neighbors of the current patch and computes their number of facial neighbors inside the
// given patch and the distance of their centers to the dipole position. Then the element with the most facial neighbors with patch elements is added.
// If multiple elements share the same number of facial neighbors in the patch, the one with the closest distance to the dipole is taken.
// I find this implementation quite unintuitive. The criterion to exclude all elements that have a facial neighbor outside gray matter seems
// somewhat arbitrary, and is contrary to the description in the paper and is contrary to how we handle the restriction in the Venant case.
// Furthermore, it seems quite arbitrary to choose a different approach for the first patch extension to the facial neighboors of the dipole
// element than for the later extensions. Furthermore, I find it unintuitive to prefer the number of facial neighbors in the set over the
// distance to the dipole. I think it might be better to keep the source configuration as focal as possible, in contrast to potentially
// adding a few more functions to the source configuration.
// Since the paper description and the implementation thus seem to diverge, I took the creative freedom to implement a custom heuristic for
// dynamic source configuration construction. The implemented heuristic works as follows.
// We first add the dipole element to create an initial patch. Then we iteratively apply the following algorithm.
// Let C be the candidate set, consisting of all elements that
//  1) share a face with the current patch,
//  2) are in the same compartment as the dipole, and
//  3) are not contained in the current patch.
// For each elements in this set, we compute the distance of the center of the element to the dipole. We then add the closest element
// in the candidate set to the patch and update the candidate set. This is repeated until the desired number of elements is reached.
template<class VolumeConductor, class GridFunctionSpace>
class DynamicSourceConfigurationBuilder
  : public SourceConfigurationBuilderInterface<TwoVerticesHdivFunction<GridFunctionSpace>>
{
public:
  using BasisFunction = TwoVerticesHdivFunction<GridFunctionSpace>;
  using DipoleType = typename BasisFunction::DipoleType;
  using ElementType = typename GridFunctionSpace::Traits::GridView::template Codim<0>::Entity;
  using GridView = typename GridFunctionSpace::Traits::GridViewType;
  using Scalar = typename GridView::ctype;
  enum {dim = GridFunctionSpace::Traits::GridViewType::dimension};
  using Index = typename GridView::IndexSet::IndexType;
  using ElementSeed = typename GridView::Grid::template Codim<0>::EntitySeed;
  using RankedElement = std::pair<ElementSeed, Scalar>;
  enum {EDGES_PER_ELEMENT = 6};
  enum {ADDITIONAL_BASIS_FUNCTIONS_PER_ELEMENT_ESTIMATE = 4};

  DynamicSourceConfigurationBuilder(std::shared_ptr<const VolumeConductor> volumeConductor,
                                    const GridFunctionSpace& gfs,
                                    const Dune::ParameterTree& config)
    : volumeConductor_(volumeConductor)
    , gfs_(gfs)
    , edgeMapper_(volumeConductor_->gridView(), Dune::mcmgLayout(Dune::Dim<1>()))
    , elementMapper_(volumeConductor_->gridView(), Dune::mcmgElementLayout())
    , n_elems_(config.get<size_t>("n_elems"))
    , basisSizeEstimate_(EDGES_PER_ELEMENT + (n_elems_ - 1) * ADDITIONAL_BASIS_FUNCTIONS_PER_ELEMENT_ESTIMATE)
    , ensureFocalFluxes_(config.get<bool>("ensureFocalFluxes", true))
  {
    if(n_elems_ < 1) {
      DUNE_THROW(Dune::Exception, "please set n_elems to a value >= 1");
    }
  }

  std::vector<BasisFunction> build(const ElementType& initialElement, const DipoleType& dipole) const override
  {
    std::vector<BasisFunction> basisFunctions;
    basisFunctions.reserve(basisSizeEstimate_);

    auto ref = referenceElement(initialElement.geometry());
    std::set<Index> usedEdges;
    std::set<Index> patchElements;

    // first add internal edges
    patchElements.insert(elementMapper_.index(initialElement));
    for(int edge_index = 0; edge_index < ref.size(dim - 1); ++edge_index) {
      auto vertexIndexRange = ref.subEntities(edge_index, dim - 1, dim);
      auto vertexIndexIterator = vertexIndexRange.begin();
      auto index_1 = *(vertexIndexIterator++);
      auto index_2 = *(vertexIndexIterator);
      basisFunctions.emplace_back(gfs_, initialElement.template subEntity<dim>(index_1), initialElement.template subEntity<dim>(index_2));
      usedEdges.insert(edgeMapper_.subIndex(initialElement, edge_index, dim - 1));
    }

    if(n_elems_ == 1) return basisFunctions;

    // initialize candidate set
    auto compareOnSecondEntry = [](const RankedElement& l, const RankedElement& r) {return l.second > r.second;};
    std::priority_queue<RankedElement, std::vector<RankedElement>, decltype(compareOnSecondEntry)> candidates(compareOnSecondEntry);

    const auto& grid = gfs_.gridView().grid();
    auto dipolePosition = dipole.position();
    auto initialConductivity = volumeConductor_->tensor(initialElement);
    Scalar maxFacialNeighborDistance = 0.0;

    for(const auto& intersection : intersections(volumeConductor_->gridView(), initialElement)) {
      if(!intersection.neighbor()) continue;

      auto outsideEntity = intersection.outside();
      if(volumeConductor_->tensor(outsideEntity) != initialConductivity) continue;

      auto distanceSquared = (outsideEntity.geometry().center() - dipolePosition).two_norm2();
      if (distanceSquared > maxFacialNeighborDistance) maxFacialNeighborDistance = distanceSquared;
      candidates.emplace(outsideEntity.seed(), distanceSquared);
    }

    // dynamically construct patch
    size_t current_size = 1;
    while(!candidates.empty()) {
      auto nextElement = grid.entity(candidates.top().first);
      candidates.pop();

      // first add face basis functions
      for(const auto& intersection : intersections(volumeConductor_->gridView(), nextElement)) {
        if(!intersection.neighbor()) continue;
        auto outsideEntity = intersection.outside();
        if(patchElements.count(elementMapper_.index(outsideEntity))) {
          auto insideVertexIndex = dim - intersection.indexInInside();
          auto outsideVertexIndex = dim - intersection.indexInOutside();
          basisFunctions.emplace_back(gfs_, nextElement.template subEntity<dim>(insideVertexIndex), outsideEntity.template subEntity<dim>(outsideVertexIndex));
        }
      }

      // now add edge basis functions
      auto nextRef = referenceElement(nextElement.geometry());
      for(int edge_index = 0; edge_index < nextRef.size(dim - 1); ++edge_index) {
        if(usedEdges.count(edgeMapper_.subIndex(nextElement, edge_index, dim - 1))) continue;
        auto vertexIndexRange = nextRef.subEntities(edge_index, dim - 1, dim);
        auto vertexIndexIterator = vertexIndexRange.begin();
        auto index_1 = *(vertexIndexIterator++);
        auto index_2 = *(vertexIndexIterator);
        basisFunctions.emplace_back(gfs_, nextElement.template subEntity<dim>(index_1), nextElement.template subEntity<dim>(index_2));
        usedEdges.insert(edgeMapper_.subIndex(nextElement, edge_index, dim - 1));
      }

      // now update candidate queue and patch
      ++current_size;
      if(current_size == n_elems_) break;

      patchElements.insert(elementMapper_.index(nextElement));
      for(const auto& intersection : intersections(volumeConductor_->gridView(), nextElement)) {
        auto outsideEntity = intersection.outside();
        if(!patchElements.count(elementMapper_.index(outsideEntity)) && volumeConductor_->tensor(outsideEntity) == initialConductivity) {
          auto distanceSquared = (outsideEntity.geometry().center() - dipolePosition).two_norm2();
          // during numerical experiments we saw that the inclusion of elements which are far removed from the source can degrade the accuracy
          // of the forward simulation. If the corresponding flag is set we thus only include elements whose center is close to the dipole
          // position, where the cutoff for "closeness" is defined by the maximal distance (squared) of the initial candidate set to the dipole position.
          if(ensureFocalFluxes_ && distanceSquared > maxFacialNeighborDistance) continue; 
          candidates.emplace(outsideEntity.seed(), distanceSquared);
        }
      }
    }

    return basisFunctions;
  }

private:
  std::shared_ptr<const VolumeConductor> volumeConductor_;
  const GridFunctionSpace& gfs_;
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> edgeMapper_;
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper_;
  const size_t n_elems_;
  const size_t basisSizeEstimate_;
  bool ensureFocalFluxes_;
};

template <class BasisFunction>
std::vector<typename BasisFunction::Scalar> PBOInterpolation(const typename BasisFunction::DipoleType& dipole, const std::vector<BasisFunction>& basisFunctions)
{
  enum {dim = BasisFunction::dim};
  using MappedArray = Eigen::Map<const Eigen::Matrix<typename BasisFunction::Scalar, dim, 1>>;
  size_t numberBasisFunctions = basisFunctions.size();
  const auto& dipolePosition = dipole.position();

  // assemble lhsMatrix = [D, Q^t; Q, 0]. We refer to Pursiainen et. al., 2016 for details
  Eigen::MatrixXd lhsMatrix(numberBasisFunctions + dim, numberBasisFunctions + dim);
  lhsMatrix.topLeftCorner(numberBasisFunctions, numberBasisFunctions) = Eigen::MatrixXd::Zero(numberBasisFunctions, numberBasisFunctions);
  lhsMatrix.bottomRightCorner<dim, dim>() = Eigen::MatrixXd::Zero(dim, dim);

  for(size_t i = 0; i < numberBasisFunctions; ++i) {
    auto syntheticDipole = basisFunctions[i].syntheticDipole();

    // D portion
    auto diff = syntheticDipole.position() - dipolePosition;
    lhsMatrix(i, i) = diff.two_norm2();

    // Q, Q^t portion
    MappedArray syntheticMomentMap(syntheticDipole.moment().data());
    lhsMatrix.col(i).tail<dim>() = syntheticMomentMap;
    lhsMatrix.row(i).tail<dim>() = syntheticMomentMap.transpose();
  }

  // assamble rhs vector
  Eigen::VectorXd rhs(numberBasisFunctions + dim);
  rhs.head(numberBasisFunctions) = Eigen::VectorXd::Zero(numberBasisFunctions);
  rhs.tail<dim>() = MappedArray(dipole.moment().data());

  Eigen::VectorXd solution = lhsMatrix.colPivHouseholderQr().solve(rhs);

  return std::vector<typename BasisFunction::Scalar>(solution.data(), solution.data() + numberBasisFunctions);
} //PBOInterpolation

template <class BasisFunction>
std::vector<typename BasisFunction::Scalar> MPOInterpolation(const typename BasisFunction::DipoleType& dipole, const std::vector<BasisFunction>& basisFunctions, const typename BasisFunction::Scalar referenceLength, bool ensureSmallFluxes)
{
  enum {dim = BasisFunction::dim};
  using MappedArray = Eigen::Map<const Eigen::Matrix<typename BasisFunction::Scalar, dim, 1>>;
  // even though values of enumerations are implicitly convertible to integral types, using them in constructors of Eigen matrices does not compile for
  // Eigen version <= 3.3.7. Once we are using a more recent Eigen version, this could also be replaced by an enum. See also
  // https://gitlab.com/libeigen/eigen/-/commit/441b3511de7d0a6930dc7643943aa5c6082ab9e1
  static constexpr size_t MPO_INTERPOLATION_CONDITIONS = dim * (dim + 1);
  size_t numberBasisFunctions = basisFunctions.size();
  const auto& dipolePosition = dipole.position();

  // assemble Matrix M. We refer to Pursiainen et. al., 2016 for details
  Eigen::MatrixXd lhsMatrix(MPO_INTERPOLATION_CONDITIONS, numberBasisFunctions);
  for(size_t j = 0; j < numberBasisFunctions; ++j) {
    auto syntheticDipole = basisFunctions[j].syntheticDipole();
    auto scaledDiff = (syntheticDipole.position() - dipolePosition) / referenceLength;
    MappedArray syntheticMomentMap(syntheticDipole.moment().data());
    lhsMatrix.col(j).head<dim>() = syntheticMomentMap;
    for(size_t k = 0; k < dim; ++k) {
      lhsMatrix.col(j).segment<dim>((k + 1) * dim) = scaledDiff[k] * syntheticMomentMap;
    }
  }

  // assemble right hand side
  Eigen::VectorXd rhs(MPO_INTERPOLATION_CONDITIONS);
  rhs.head<dim>() = MappedArray(dipole.moment().data());
  rhs.tail<dim * dim>() = Eigen::VectorXd::Zero(dim * dim);

  // compute interpolation coefficients as the minimum norm solution of lhsMatrix x = rhs
  // NOTE : Passing the computation flags as function arguments is deprecated in Eigen 3.4.
  // But since at the point of programming this function the Ubuntu package repository installs Eigen 3.3,
  // which does not support passing these arguments as template parameters, it will stay for now.
  Eigen::VectorXd solution = lhsMatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);

  // in our multilayer sphere experiments the lhs matrix seemed to always be rank deficient. The SVD then sometimes included very
  // small singular values in the pseudoinverse, leading to what one would call "ghost sources" in the inverse problem. By this we mean
  // very large currents which cancel out to produce small currents. Imagine writing (0, 1) = 100 * (1, 0.01) - 100 * (1, 0).
  // especially if the currents are close to a conductivity jump this can massively degrade the forward modelling accuracy, and we thus want
  // to exclude these interpolations. To this end we want to only allow interpolations where the dipole moment is additively build up from
  // small fluxes, instead of arising by the cancellation of large fluxes. Our approach is to label fluxes as "small" as long as their synthetic
  // dipole moment is not larger in norm than that of the dipole they aim to interpolate. This is quite strict and excludes some "valid" interpolations
  // producing accurate results, but it also prevents certain cases where the interpolations completely fails. We are thus trading potential accuracy 
  // for stability.
  if(ensureSmallFluxes && solution.array().abs().maxCoeff() > dipole.moment().two_norm()) return PBOInterpolation(dipole, basisFunctions);

  return std::vector<typename BasisFunction::Scalar>(solution.data(), solution.data() + numberBasisFunctions);
} //MPOInterpolation

/*
 * Assembling the FEM right hand side using the Whitney source model consists of three steps.
 *
 * 1): Decide on a set of H(div) functions w_1, ..., w_k. We call such a selection a
 *     source configuration
 *
 * 2): Interpolate a given dipole M \delta_{x_0} into the space spanned by w_1, ..., w_k, which
 *     means computing coefficients c_1, ..., c_k so that
 *     M \delta_{x_0} \approx \sum_{l = 1}^k c_l w_l
 *
 * 3): Assemble the FEM right hand side corresponding to the function \sum_{l = 1}^k c_l w_l
 *
 */
template<class VolumeConductor, class GridFunctionSpace, class Vector>
class WhitneySourceModel
  : public SourceModelBase<typename GridFunctionSpace::Traits::GridViewType, Vector> {
public:
  using BaseType = SourceModelBase<typename GridFunctionSpace::Traits::GridViewType, Vector>;
  using DipoleType = typename BaseType::DipoleType;
  using GridView = typename GridFunctionSpace::Traits::GridViewType;
  enum {dim = GridView::dimension};
  using Scalar = typename GridView::ctype;
  using SearchType = typename BaseType::SearchType;
  using DOFVector = typename BaseType::VectorType;
  enum class Interpolation {PBO, MPO, AUTO};
  using BasisFunction = TwoVerticesHdivFunction<GridFunctionSpace>;

  WhitneySourceModel(std::shared_ptr<const VolumeConductor> volumeConductor,
                     const GridFunctionSpace& gfs,
                     std::shared_ptr<const SearchType> search,
                     const Dune::ParameterTree& config)
    : BaseType(search)
    , sourceConfigBuilder_(createSourceConfigurationBuilder(volumeConductor, gfs, config))
    , interpolation_(interpolationFromString(config.get<std::string>("interpolation")))
    , referenceLength_(interpolation_ == Interpolation::PBO ? 0.0 : config.get<Scalar>("referenceLength"))
    , ensureSmallFluxes_(interpolation_ == Interpolation::PBO ? false : config.get<bool>("ensureSmallFluxes", true))
  {
  }

  void assembleRightHandSide(DOFVector& dofVector) const override
  {
    // assemble source configuration
    auto basisFunctions = sourceConfigBuilder_->build(this->dipoleElement(), this->dipole());

    // interpolate dipole into source configuration
    std::vector<Scalar> coefficients;
    if(interpolation_ == Interpolation::PBO) {
      coefficients = PBOInterpolation(this->dipole(), basisFunctions);
    }
    else if(interpolation_ == Interpolation::MPO) {
      coefficients = MPOInterpolation(this->dipole(), basisFunctions, referenceLength_, ensureSmallFluxes_);
    }
    else {
      // In the MPO approach, the number of constraints is <= the number of degrees of freedoms iff the condition below is fulfilled. In Pursiainen et. al., 2016,
      // it is noted that in this not overdetermined setting the MPO interpolation approach seems favourable. We follow this suggestion.
      coefficients = basisFunctions.size() >= dim * (dim + 1) ? MPOInterpolation(this->dipole(), basisFunctions, referenceLength_, ensureSmallFluxes_)
                                                              : PBOInterpolation(this->dipole(), basisFunctions);
    }

    // assemble right hand side of interpolated dipole
    for(size_t i = 0; i < basisFunctions.size(); ++i) {
      basisFunctions[i].accumulate(dofVector, coefficients[i]);
    }
  }

private:
  const std::unique_ptr<const SourceConfigurationBuilderInterface<BasisFunction>> sourceConfigBuilder_;
  const Interpolation interpolation_;
  const Scalar referenceLength_;
  const bool ensureSmallFluxes_;

  Interpolation interpolationFromString(const std::string& value) {
    if(value == "PBO") {
      return Interpolation::PBO;
    }
    else if(value == "MPO") {
      return Interpolation::MPO;
    }
    else if(value == "auto") {
      return Interpolation::AUTO;
    }
    else {
      DUNE_THROW(Dune::Exception, "interpolation has invalid value " << value << ", please choose one from the list {PBO, MPO, auto}");
    }
  }

  std::unique_ptr<SourceConfigurationBuilderInterface<BasisFunction>> createSourceConfigurationBuilder(std::shared_ptr<const VolumeConductor> volumeConductor,
                                                                                                       const GridFunctionSpace& gfs,
                                                                                                       const Dune::ParameterTree& config)
  {
    std::string value = config.get<std::string>("sourceConfiguration");
    if(value == "fixed") {
      return std::make_unique<FixedSourceConfigurationBuilder<VolumeConductor, GridFunctionSpace>>(volumeConductor, gfs, config);
    }
    else if(value == "dynamic") {
      return std::make_unique<DynamicSourceConfigurationBuilder<VolumeConductor, GridFunctionSpace>>(volumeConductor, gfs, config);
    }
    else {
      DUNE_THROW(Dune::Exception, "sourceConfiguration has invalid value " << value <<", please choose one from the list {fixed, dynamic}");
    }
  }
}; // class WhitneySourceModel

} // namespace duneuro

#endif // DUNEURO_WHITNEY_SOURCE_MODEL_HH
