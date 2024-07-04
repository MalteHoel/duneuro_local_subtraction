#ifndef DUNEURO_SOURCE_SPACE_FACTORY_HH
#define DUNEURO_SOURCE_SPACE_FACTORY_HH

#include <memory>

#include <dune/common/parametertree.hh>

#if HAVE_DUNE_UDG
#include <dune/udg/pdelab/subtriangulation.hh>
#endif

#include <duneuro/common/exceptions.hh>
#include <duneuro/driver/feature_manager.hh>

namespace duneuro
{
  struct sourceSpaceFactory {
#if HAVE_DUNE_UDG
// basic source space for unfitted grid. Iterates over fundamental mesh with at least "thresh"% gm proportion, places source in center
    template< typename GFS, typename ST>
    std::vector<std::vector<double>>                                                
    createUnfitted(const GFS& gfs, const ST& subTriangulation, const Dune::ParameterTree& config)
    {
      using GV = typename GFS::Traits::GridViewType;
      using UST = Dune::PDELab::UnfittedSubTriangulation<GV>;

      auto gridView = gfs.gridView();
      auto compInd = config.get<int>("sourceCompartmentIndex");
     // auto gridRes = config.get<int>("sourceGridResolution");
      std::vector<std::vector<double>> outputGrid;
      auto thresh = 0.2;
      for (const auto& element : Dune::elements(gridView)) { 
        // skip non-gm elements
        if (!subTriangulation.isHostCell(element, compInd)) {
          continue; }
        auto vol = element.geometry().volume(); 
        auto gmVol = 0;
        UST ust(subTriangulation.gridView(), subTriangulation);
        ust.create(element);
        for (const auto& ep : ust) {
          // skip snippets outside gm
          if (ep.domainIndex() != compInd){
            continue; 
          }
          gmVol += ep.geometry().volume();
          // check if threshold condition for gm proportion is met yet
          if (gmVol/vol < thresh){
            continue;
          }
          // add center of latest snippet with x,y,z dipole moment to grid
          std::vector<double> source(6);
          for (int i = 0 ; i< ep.geometry().center().size(); i++){
            source[i] = ep.geometry().center()[i];
          }
          source[3] = 1;
          source[4] = 0;
          source[5] = 0;
          outputGrid.push_back(source);
          source[3] = 0;
          source[4] = 1;
          outputGrid.push_back(source);
          source[4] = 0;
          source[5] = 1;
          outputGrid.push_back(source);
          break;
        }
      }
      return outputGrid;
    }
#endif

    std::vector<std::vector<double>> 
    createFitted(const Dune::ParameterTree& config)
    {
    DUNE_THROW(Dune::NotImplemented, "currently not implemented");
    }
  };
}

#endif // DUNEURO_SOURCE_SPACE_FACTORY_HH
