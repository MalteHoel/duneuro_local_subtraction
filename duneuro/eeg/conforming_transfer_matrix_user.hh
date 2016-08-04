#ifndef DUNEURO_CONFORMING_TRANSFER_MATRIX_USER_HH
#define DUNEURO_CONFORMING_TRANSFER_MATRIX_USER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dg_solver.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/sparse_vector_container.hh>
#include <duneuro/common/vector_density.hh>
#include <duneuro/eeg/dg_source_model_factory.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class S>
  struct ConformingTransferMatrixUserTraits {
    static const unsigned int dimension = S::Traits::dimension;
    using VolumeConductor = typename S::Traits::VolumeConductor;
    using EEGForwardSolver = S;
    using DenseRHSVector = typename EEGForwardSolver::Traits::RangeDOFVector;
    using SparseRHSVector = SparseVectorContainer<typename DenseRHSVector::ContainerIndex,
                                                  typename DenseRHSVector::ElementType>;
    using CoordinateField = typename VolumeConductor::ctype;
    using Coordinate = Dune::FieldVector<CoordinateField, dimension>;
    using DipoleType = Dipole<CoordinateField, dimension>;
    using DomainField = typename EEGForwardSolver::Traits::DomainDOFVector::field_type;
  };

  template <class S, class SMF>
  class ConformingTransferMatrixUser
  {
  public:
    using Traits = ConformingTransferMatrixUserTraits<S>;

    ConformingTransferMatrixUser(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                                 std::shared_ptr<typename Traits::EEGForwardSolver> solver)
        : volumeConductor_(volumeConductor), solver_(solver)
    {
    }

    ConformingTransferMatrixUser(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                                 const Dune::ParameterTree& config)
        : ConformingTransferMatrixUser(
              volumeConductor,
              std::make_shared<typename Traits::EEGForwardSolver>(volumeConductor, config), config)
    {
    }

    void postProcessPotential(const typename Traits::DipoleType& dipole,
                              const std::vector<typename Traits::Coordinate>& projectedElectrodes,
                              std::vector<typename Traits::DomainField>& potential,
                              const Dune::ParameterTree& config)
    {
      auto density = source_model_default_density(config.sub("source_model"));
      if (density == VectorDensity::sparse) {
        auto sourceModel = SMF::template createSparse<typename Traits::SparseRHSVector>(
            volumeConductor_, *solver_, config.sub("source_model"));
        sourceModel->postProcessSolution(dipole, projectedElectrodes, potential);
      } else {
        auto sourceModel = SMF::template createDense<typename Traits::DenseRHSVector>(
            volumeConductor_, *solver_, config.sub("source_model"));
        sourceModel->postProcessSolution(dipole, projectedElectrodes, potential);
      }
    }

    template <class M>
    std::vector<typename Traits::DomainField>
    solve(const M& transferMatrix, const typename Traits::DipoleType& dipole,
          const Dune::ParameterTree& config, DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      std::vector<typename Traits::DomainField> result;
      auto density = source_model_default_density(config.sub("source_model"));
      if (density == VectorDensity::sparse) {
        dataTree.set("density", "sparse");
        result = solveSparse(transferMatrix, dipole, config);
      } else {
        dataTree.set("density", "dense");
        result = solveDense(transferMatrix, dipole, config);
      }
      dataTree.set("time", timer.elapsed());
      return result;
    }

    template <class M>
    std::vector<typename Traits::DomainField> solveSparse(const M& transferMatrix,
                                                          const typename Traits::DipoleType& dipole,
                                                          const Dune::ParameterTree& config) const
    {
      using SVC = typename Traits::SparseRHSVector;
      SVC rhs;
      auto sourceModel = SMF::template createSparse<typename Traits::SparseRHSVector>(
          volumeConductor_, *solver_, config.sub("source_model"));
      sourceModel->assembleRightHandSide(dipole, rhs);

      const auto blockSize =
          Traits::EEGForwardSolver::Traits::FunctionSpace::GFS::Traits::Backend::blockSize;

      if (blockSize == 1) {
        return matrix_sparse_vector_product(transferMatrix, rhs,
                                            [](const typename SVC::Index& c) { return c[0]; });
      } else {
        return matrix_sparse_vector_product(
            transferMatrix, rhs,
            [blockSize](const typename SVC::Index& c) { return c[1] * blockSize + c[0]; });
      }
    }

    template <class M>
    std::vector<typename Traits::DomainField> solveDense(const M& transferMatrix,
                                                         const typename Traits::DipoleType& dipole,
                                                         const Dune::ParameterTree& config) const
    {
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(*solver_, 0.0);
      } else {
        *denseRHSVector_ = 0.0;
      }
      auto sourceModel = SMF::template createDense<typename Traits::DenseRHSVector>(
          volumeConductor_, *solver_, config.sub("source_model"));
      sourceModel->assembleRightHandSide(dipole, *denseRHSVector_);

      return matrix_dense_vector_product(transferMatrix,
                                         Dune::PDELab::Backend::native(*denseRHSVector_));
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::EEGForwardSolver> solver_;
    mutable std::shared_ptr<typename Traits::DenseRHSVector> denseRHSVector_;
  };
}

#endif // DUNEURO_CONFORMING_TRANSFER_MATRIX_USER_HH
