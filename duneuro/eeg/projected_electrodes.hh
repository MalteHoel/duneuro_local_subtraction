#ifndef DUNEURO_PROJECTEDELECTRODES_HH
#define DUNEURO_PROJECTEDELECTRODES_HH

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <ostream>
#include <vector>

#if HAVE_DUNE_UDG
#include <dune/udg/pdelab/multiphaseoperator.hh>
#include <dune/udg/pdelab/operator.hh>
#include <dune/udg/pdelab/subtriangulation.hh>
#include <duneuro/common/kdtree.hh>
#endif
//#include <dune/biomag/localoperator/boundaryprojection.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class GV>
  class ProjectedElectrodes : public ElectrodeProjectionInterface<GV>
  {
    typedef typename GV::template Codim<0>::Entity Element;

  public:
    typedef typename GV::ctype ctype;
    enum { dim = GV::dimension };

    using Projection = ProjectedElectrode<GV>;

    ProjectedElectrodes(const std::vector<Dune::FieldVector<ctype, dim>>& electrodes, const GV& gv,
                        DataTree dataTree = DataTree())
        : gridView_(gv)
    {
      Dune::Timer timer;
      std::vector<std::pair<Projection, ctype>> minDistance;
      for (std::size_t i = 0; i < electrodes.size(); ++i) {
        minDistance.push_back(std::make_pair(
            Projection{gridView_.template begin<0>(), Dune::FieldVector<ctype, dim>(0.0)},
            std::numeric_limits<ctype>::max()));
      }
      std::cout << "\n";
      for (const auto& element : elements(gridView_)) {
        if (!element.hasBoundaryIntersections())
          continue;
        const auto& eg = element.geometry();
        const auto& elementReference = Dune::ReferenceElements<ctype, dim>::general(eg.type());
        for (unsigned int i = 0; i < electrodes.size(); ++i) {
          auto local = eg.local(electrodes[i]);
          for (unsigned int i = 0; i < dim; ++i) {
            local[i] = std::max(std::min(local[i], 1.0), 0.0);
          }
          // if (elementReference.checkInside(local)) {
          auto diff = electrodes[i];
          diff -= eg.global(local);
          auto diff2n = diff.two_norm();
          if (diff2n < minDistance[i].second) {
            minDistance[i].first = Projection{element, local};
            minDistance[i].second = diff2n;
          }
          //}
        }
      }
      ctype maxdiff = 0.0;
      for (unsigned int i = 0; i < electrodes.size(); ++i) {
        maxdiff = std::max(maxdiff, minDistance[i].second);
        if (minDistance[i].second == std::numeric_limits<ctype>::max()) {
          DUNE_THROW(Dune::Exception, "electrode " << i << " at " << electrodes[i]
                                                   << " could not be projected");
        }
        projections_.push_back(minDistance[i].first);
      }
      dataTree.set("max_distance", maxdiff);
      dataTree.set("time", timer.elapsed());
    }

#if HAVE_DUNE_UDG
    template <class GFS, class ST>
    ProjectedElectrodes(const std::vector<Dune::FieldVector<ctype, dim>>& electrodes,
                        const GFS& gfs, const ST& subTriangulation, DataTree dataTree = DataTree())
        : gridView_(gfs.gridView())
    {
      Dune::Timer timer;
      KDTreeElementSearch<GV> search(gridView_);
      std::vector<int> insideHostCell(electrodes.size(),std::numeric_limits<ctype>::max());
      projections_.resize(electrodes.size());
      int elecOffset = 0;
      // use kd-search to find element that contains the electrode
      for (const auto& electrode : electrodes) {
        auto search_result = search.findEntity(electrode);
        if(!search_result.has_value()) {
          DUNE_THROW(Dune::Exception, "coordinate is outside of the grid, or grid is not convex");
        }
        const auto& element = search_result.value();
        if (subTriangulation.isHostCell(element)) {
            projections_[elecOffset] = Projection{element, element.geometry().local(electrode)};
            insideHostCell[elecOffset] = 1;
        }
        elecOffset +=1;
      }
      // for the elements that are not host to a domain, perform a grid wide search for the closest vertex,
      // identical to closest_subentity_center_electrode_projection

      std::vector<ctype> distances(projections_.size(), std::numeric_limits<ctype>::max());
      int codim = 3;

      for (const auto& element : Dune::elements(gridView_)) { 
        if (!subTriangulation.isHostCell(element)) {
          continue;}
        // skip electrodes that were already projected
        for (int elecInd = 0 ; elecInd < electrodes.size(); elecInd++){
          if (insideHostCell[elecInd] == 1) {
          continue; }
          const auto& geo = element.geometry();
          const auto& ref = Dune::ReferenceElements<ctype, GV::dimension>::general(geo.type());
          for (int i = 0; i < ref.size(codim); ++i) {
            auto local = ref.position(i, codim);
            auto diff = geo.global(local);
            diff -= electrodes[elecInd];
            auto distance = diff.two_norm();
            if (distance < distances[elecInd]) {
              distances[elecInd] = distance;
              projections_[elecInd] = Projection{element, local};
            }
          }
         // std::cout<< distances[elecInd] << std::endl;
        }

      }
         // DUNE_THROW(Dune::Exception, "element of electrode at "
         //                                 << electrode << " is not a host cell for any domain");
      dataTree.set("time", timer.elapsed());
    }
#endif

    virtual void
    setElectrodes(const std::vector<Dune::FieldVector<ctype, dim>>& electrodes) override
    {
      DUNE_THROW(Dune::Exception, "should not be called");
    }

    template <class DGF, class OutputIterator>
    void evaluateAtProjections(const DGF& dgf, OutputIterator out) const
    {
      for (std::size_t i = 0; i < projections_.size(); ++i) {
        typename DGF::Traits::RangeType y(0.0);
        dgf.evaluate(projections_[i].element, projections_[i].localPosition, y);
        *out++ = y;
      }
    }

    template <class GFS, class X, class OutputIterator>
    void evaluateAtProjections(const GFS& gfs, const X& x, OutputIterator out) const
    {
      using DGF = Dune::PDELab::DiscreteGridFunction<GFS, X>;
      evaluateAtProjections(DGF(gfs, x), out);
    }

    template <class GFS, class X>
    std::vector<typename X::field_type> evaluate(const GFS& gfs, const X& x) const
    {
      std::vector<typename X::field_type> out;
      out.reserve(projections_.size());
      evaluateAtProjections(gfs, x, std::back_inserter(out));
      std::cout << "got " << out.size() << " out"
                << "\n";
      return out;
    }

    Dune::FieldVector<ctype, dim> projection(std::size_t i) const
    {
      return element(i).geometry().global(projections_[i].localPosition);
    }

    const Element& element(std::size_t i) const
    {
      return projections_[i].element;
    }

    virtual const Projection& getProjection(std::size_t i) const override
    {
      return projections_[i];
    }

    virtual std::size_t size() const override
    {
      return projections_.size();
    }

  private:
    GV gridView_;
    std::vector<Projection> projections_;
  };
}

#endif // DUNEURO_PROJECTEDELECTRODES_HH

