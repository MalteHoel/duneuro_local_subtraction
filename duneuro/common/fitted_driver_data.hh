#ifndef DUNEURO_FITTED_DRIVER_DATA_HH
#define DUNEURO_FITTED_DRIVER_DATA_HH

#include <memory>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <duneuro/common/simple_structured_grid.hh>

namespace duneuro
{
  /**
   * \brief data storage for fitted methods
   *
   * This class stores data used for initialization of fitted methods, e.g. CG or DG. Its main use
   * is to provide mesh data which has not been read from a file but provide from outside. This is
   * used when creating a mesh from python or Matlab data.
   */
  template <int dim>
  struct FittedDriverData {
    using Coordinate = Dune::FieldVector<double, dim>;
    using Tensor = Dune::FieldMatrix<double,dim,dim>;
    //! \brief nodes of the mesh
    std::vector<Coordinate> nodes;
    //! \brief elements of the mesh, defined by the indices of their nodes in the `nodes` vector
    std::vector<std::vector<unsigned int>> elements;
    //! \brief on label for each element describing its conductivity
    std::vector<int> labels;
    //! \brief list of conductivity values, e.g. 4 entries for a 4 compartment isotropic model
    std::vector<double> conductivities;
    //! \brief list of conductivity tensors for each element
    std::vector<Tensor> tensors;
  };
}

#endif // DUNEURO_FITTED_DRIVER_DATA_HH
