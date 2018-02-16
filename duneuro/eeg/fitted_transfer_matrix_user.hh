#ifndef DUNEURO_FITTED_TRANSFER_MATRIX_USER_HH
#define DUNEURO_FITTED_TRANSFER_MATRIX_USER_HH

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
  struct FittedTransferMatrixUserTraits {
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
    using ElementSearch = KDTreeElementSearch<typename VolumeConductor::GridView>;
  };

  template <class S, class SMF>
  class FittedTransferMatrixUser
  {
  public:
    using Traits = FittedTransferMatrixUserTraits<S>;

    FittedTransferMatrixUser(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                             std::shared_ptr<typename Traits::ElementSearch> search,
                             std::shared_ptr<typename Traits::EEGForwardSolver> solver)
        : volumeConductor_(volumeConductor), search_(search), solver_(solver)
    {
    }

    void setSourceModel(const Dune::ParameterTree& config, const Dune::ParameterTree& solverConfig,
                        DataTree dataTree = DataTree())
    {
      sparseSourceModel_.reset();
      denseSourceModel_.reset();
      density_ = source_model_default_density(config);
      if (density_ == VectorDensity::sparse) {
        sparseSourceModel_ = SMF::template createSparse<typename Traits::SparseRHSVector>(
            volumeConductor_, *solver_, search_, config, solverConfig);
      } else {
        denseSourceModel_ = SMF::template createDense<typename Traits::DenseRHSVector>(
            volumeConductor_, *solver_, search_, config, solverConfig);
      }
    }

    void bind(const typename Traits::DipoleType& dipole, DataTree dataTree = DataTree())
    {
      if (density_ == VectorDensity::sparse) {
        if (!sparseSourceModel_) {
          DUNE_THROW(Dune::Exception, "source model not set");
        }
        sparseSourceModel_->bind(dipole, dataTree);
      } else {
        if (!denseSourceModel_) {
          DUNE_THROW(Dune::Exception, "source model not set");
        }
        denseSourceModel_->bind(dipole, dataTree);
      }
    }

    void postProcessPotential(const std::vector<typename Traits::Coordinate>& projectedElectrodes,
                              std::vector<typename Traits::DomainField>& potential)
    {
      if (projectedElectrodes.size() != potential.size()) {
        DUNE_THROW(duneuro::IllegalArgumentException,
                   "number of electrodes ("
                       << projectedElectrodes.size()
                       << ") does not match number of entries in the potential vector ("
                       << potential.size() << ")");
      }
      if (density_ == VectorDensity::sparse) {
        if (!sparseSourceModel_) {
          DUNE_THROW(Dune::Exception, "source model not set");
        }
        sparseSourceModel_->postProcessSolution(projectedElectrodes, potential);
      } else {
        if (!denseSourceModel_) {
          DUNE_THROW(Dune::Exception, "source model not set");
        }
        denseSourceModel_->postProcessSolution(projectedElectrodes, potential);
      }
    }

    template <class M>
    std::vector<typename Traits::DomainField> solve(const M& transferMatrix,
                                                    DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      std::vector<typename Traits::DomainField> result;
      if (density_ == VectorDensity::sparse) {
        dataTree.set("density", "sparse");
        result = solveSparse(transferMatrix);
      } else {
        dataTree.set("density", "dense");
        result = solveDense(transferMatrix);
      }
      dataTree.set("time", timer.elapsed());
      return result;
    }

    template <class M>
    std::vector<typename Traits::DomainField> solveSparse(const M& transferMatrix) const
    {
      using SVC = typename Traits::SparseRHSVector;
      SVC rhs;
      if (!sparseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
      sparseSourceModel_->assembleRightHandSide(rhs);

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
    std::vector<typename Traits::DomainField> solveDense(const M& transferMatrix) const
    {
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(*solver_, 0.0);
      } else {
        *denseRHSVector_ = 0.0;
      }
      if (!denseSourceModel_) {
        DUNE_THROW(Dune::Exception, "source model not set");
      }
      denseSourceModel_->assembleRightHandSide(*denseRHSVector_);

      return matrix_dense_vector_product(transferMatrix,
                                         Dune::PDELab::Backend::native(*denseRHSVector_));
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::ElementSearch> search_;
    std::shared_ptr<typename Traits::EEGForwardSolver> solver_;
    VectorDensity density_;
    std::shared_ptr<SourceModelInterface<typename Traits::DomainField, Traits::dimension,
                                         typename Traits::SparseRHSVector>>
        sparseSourceModel_;
    std::shared_ptr<SourceModelInterface<typename Traits::DomainField, Traits::dimension,
                                         typename Traits::DenseRHSVector>>
        denseSourceModel_;
    mutable std::shared_ptr<typename Traits::DenseRHSVector> denseRHSVector_;
  };
}

#endif // DUNEURO_FITTED_TRANSFER_MATRIX_USER_HH
