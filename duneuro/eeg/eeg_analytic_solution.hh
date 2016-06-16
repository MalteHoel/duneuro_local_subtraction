#ifndef DUNEURO_ANALYTICSOLUTION_HH
#define DUNEURO_ANALYTICSOLUTION_HH

#include <cassert>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/legacy/analytic_solution.hh>

namespace duneuro
{
  /** compute the analytic solution of the eeg-forward problem in a multi
   *  layer sphere model
   */
  template <class ctype>
  class AnalyticSolution
  {
  public:
    enum { dim = 3 };

    AnalyticSolution(std::vector<ctype> sphereRadii, Dune::FieldVector<ctype, dim> sphereCenter,
                     std::vector<ctype> conductivities,
                     std::vector<Dune::FieldVector<ctype, dim>> electrodePositions)
        : sphereRadii_(sphereRadii)
        , sphereCenter_(sphereCenter)
        , conductivities_(conductivities)
        , electrodePositions_(electrodePositions)
    {
      if (sphereRadii_.size() != conductivities_.size()) {
        DUNE_THROW(Dune::Exception, "number of spheres does not match number of conductivities ("
                                        << sphereRadii_.size() << " vs " << conductivities_.size()
                                        << ")");
      }
    }

    AnalyticSolution(const Dune::ParameterTree& config,
                     std::vector<Dune::FieldVector<ctype, dim>> electrodePositions)
        : sphereRadii_(config.get<std::vector<ctype>>("radii"))
        , sphereCenter_(config.get<Dune::FieldVector<ctype, dim>>("center"))
        , conductivities_(config.get<std::vector<ctype>>("conductivities"))
        , electrodePositions_(electrodePositions)
    {
      if (sphereRadii_.size() != conductivities_.size()) {
        DUNE_THROW(Dune::Exception, "number of spheres does not match number of conductivities ("
                                        << sphereRadii_.size() << " vs " << conductivities_.size()
                                        << ")");
      }
    }

    std::vector<ctype> operator()(const Dipole<ctype, dim>& dipole) const
    {
      return Legacy::analytic_solution(sphereRadii_.size(), sphereRadii_, sphereCenter_,
                                       conductivities_, electrodePositions_, dipole.moment(),
                                       dipole.position());
    }

  private:
    unsigned int numberOfLayers_;
    std::vector<ctype> sphereRadii_;
    Dune::FieldVector<ctype, dim> sphereCenter_;
    std::vector<ctype> conductivities_;
    std::vector<Dune::FieldVector<ctype, dim>> electrodePositions_;
  };

  template <class ctype, int dim>
  Dune::DynamicMatrix<ctype>
  compute_analytic_solution(const std::vector<Dune::FieldVector<ctype, dim>>& electrodes,
                            const std::vector<Dipole<ctype, dim>>& dipoles,
                            const Dune::ParameterTree& config, DataTree output = DataTree())
  {
    Dune::Timer timer;
    Dune::DynamicMatrix<ctype> leadFieldMatrix(electrodes.size(), dipoles.size());

    AnalyticSolution<ctype> as(config, electrodes);

    std::vector<std::vector<ctype>> solutionsAtElectrodes;
    std::transform(dipoles.begin(), dipoles.end(), std::back_inserter(solutionsAtElectrodes), as);

    for (unsigned int i = 0; i < solutionsAtElectrodes.size(); ++i)
      for (unsigned int j = 0; j < solutionsAtElectrodes[i].size(); ++j)
        leadFieldMatrix[j][i] = solutionsAtElectrodes[i][j];

    duneuro::subtract_column_means(*leadFieldMatrix);
    output.set("time", timer.elapsed());
    return leadFieldMatrix;
  }

  template <class ctype, int dim>
  std::vector<ctype>
  compute_analytic_solution(const std::vector<Dune::FieldVector<ctype, dim>>& electrodes,
                            const Dipole<ctype, dim>& dipole, const Dune::ParameterTree& config,
                            DataTree output = DataTree())
  {
    return AnalyticSolution<ctype>(config, electrodes)(dipole);
  }
}

#endif // DUNEURO_ANALYTICSOLUTION_HH
