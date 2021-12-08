#ifndef DUNEURO_ELEMENT_PATCH_ASSEMBLER_HH
#define DUNEURO_ELEMENT_PATCH_ASSEMBLER_HH

#include <memory>
#include <vector>

#include <dune/common/parametertree.hh>
#include <duneuro/common/element_patch.hh>

namespace duneuro
{
  /**
     \brief Assembler for functionals with local support on a patch
   */
  template <typename VC, typename FS> // , FittedSolverType solverType>
  class ElementPatchAssembler
  {
    using GV = typename VC::GridView;
  public:
    using VolumenConductor = VC;
    using GridView = GV;
    using Coordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using Element = typename GV::template Codim<0>::Entity;
    using Intersection = typename GV::Intersection;
    using SearchType = KDTreeElementSearch<GV>;
    using LFS = Dune::PDELab::LocalFunctionSpace<typename FS::GFS>;
    using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;

    ElementPatchAssembler(std::shared_ptr<const VC> volumeConductor,
                          std::shared_ptr<const FS> fs,
                          std::shared_ptr<const SearchType> search,
                          const Dune::ParameterTree& config)
        : search_(search)
        , volumeConductor_(volumeConductor)
        , functionSpace_(fs)
        , elementNeighborhoodMap_(std::make_shared<ElementNeighborhoodMap<typename VC::GridView>>(
              volumeConductor_->gridView()))
        , config_(config)
        , lfs_inside(functionSpace_->getGFS())
        , cache_inside(lfs_inside)
        , lfs_outside(functionSpace_->getGFS())
        , cache_outside(lfs_outside)
    {}

    void bind(const Coordinate & pos,
              DataTree dataTree = DataTree())
    {
      // create patch, according to config
      auto elementPatch = make_element_patch(volumeConductor_, elementNeighborhoodMap_,
                                             *search_, pos, config_);
      // extract patch elements
      patchElements_ = elementPatch->elements();
      dataTree.set("elements", patchElements_.size());

      // extract patch boundary intersection
      patchBoundaryIntersections_ = elementPatch->extractBoundaryIntersections();
      dataTree.set("intersections", patchBoundaryIntersections_.size());
      
      // if we use a conforming discretization, we have to extend the patch by another layer
      transitionElements_ = elementPatch->transitionElements();
      dataTree.set("transitionElements", transitionElements_.size());
    }

    template<typename Vector, typename LOP>
    void assemblePatchVolume(Vector& vector, const LOP& lop) const
    {
      assembleElementSetVolume(vector, patchElements_,
        [&lop](auto& eg, auto& lfs_inside, auto& view_inside){
          lop.lambda_patch_volume(eg, lfs_inside, view_inside);
        });
    }

    template<typename Vector, typename LOP>
    void assemblePatchBoundary(Vector& vector, const LOP& lop) const
    {
      for (const auto& is : patchBoundaryIntersections_) {
        // retrieve and bind inside
        lfs_inside.bind(is.inside());
        cache_inside.update();

        // retrieve and bind outside
        lfs_outside.bind(is.outside());
        cache_outside.update();

        // resize local vectors
        v_inside.assign(cache_inside.size(), 0.0);
        v_outside.assign(cache_outside.size(), 0.0);

        // create geometry wrapper
        Dune::PDELab::IntersectionGeometry<typename GridView::Intersection> ig(is, 0);

        // call local operator
        auto view_inside = v_inside.weightedAccumulationView(1.0);
        auto view_outside = v_outside.weightedAccumulationView(1.0);
        lop.lambda_patch_boundary(ig, lfs_inside, lfs_outside, view_inside, view_outside);

        // copy back to main vector
        for (unsigned int i = 0; i < cache_inside.size(); i++) {
          auto index = cache_inside.containerIndex(i);
          vector[index] += v_inside(lfs_inside,i);
        }
        for (unsigned int i = 0; i < cache_outside.size(); i++) {
          auto index = cache_outside.containerIndex(i);
          vector[index] += v_outside(lfs_outside,i);
        }
      }
    }
    
    template<typename Vector, typename LOP>
    void assembleTransitionVolume(Vector& vector, const LOP& lop) const
    {
      assembleElementSetVolume(vector, transitionElements_,
        [&lop](auto& eg, auto& lfs_inside, auto& view_inside){
          lop.lambda_transition_volume(eg, lfs_inside, view_inside);
        });
    }
    
    const std::vector<Element>& patchElements() {
      return patchElements_;
    }
    
    const std::vector<Intersection>& intersections() {
      return patchBoundaryIntersections_;
    }
    
    const std::vector<Element>& transitionElements() {
      return transitionElements_;
    }
    
  private:
    using FESwitch =
      Dune::FiniteElementInterfaceSwitch<typename LFS::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<typename FESwitch::Basis>;
    using RF = typename BasisSwitch::RangeField;

    mutable Dune::PDELab::LocalVector<RF> v_inside;
    mutable Dune::PDELab::LocalVector<RF> v_outside;

    template<typename Vector, typename Caller>
    void assembleElementSetVolume(Vector& vector,
      const std::vector<Element>& elements,
      Caller&& caller) const
    {
      for(const auto& element : elements) {
        // retrieve and bind inside
        lfs_inside.bind(element);
        cache_inside.update();

        // resize local vector
        v_inside.assign(cache_inside.size(), 0.0);

        // create geometry wrapper
        Dune::PDELab::ElementGeometry<Element> eg(element);

        // call local operator
        auto view_inside = v_inside.weightedAccumulationView(1.0);
        caller(eg, lfs_inside, view_inside);

        // copy back to main vector
        for (unsigned int i = 0; i < cache_inside.size(); i++) {
          auto index = cache_inside.containerIndex(i);
          vector[index] += v_inside(lfs_inside,i);
        }
      }
    }
    
    std::shared_ptr<const VC> volumeConductor_;
    std::shared_ptr<const SearchType> search_;
    std::shared_ptr<const FS> functionSpace_;
    std::shared_ptr<ElementNeighborhoodMap<typename VC::GridView>> elementNeighborhoodMap_;
    Dune::ParameterTree config_;

    mutable LFS lfs_inside;
    mutable LFSCache cache_inside;
    mutable LFS lfs_outside;
    mutable LFSCache cache_outside;
    
    std::vector<Element> patchElements_;
    std::vector<Intersection> patchBoundaryIntersections_;
    std::vector<Element> transitionElements_;
  };

}

#endif // DUNEURO_ELEMENT_PATCH_ASSEMBLER_HH
