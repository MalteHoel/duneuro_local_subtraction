#ifndef DUNEURO_TRANSFER_MATRIX_USER_HH
#define DUNEURO_TRANSFER_MATRIX_USER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/sparse_vector_container.hh>
#include <duneuro/common/vector_density.hh>
#include <duneuro/io/data_tree.hh>
#include <duneuro/eeg/electrode_projection_interface.hh>

namespace duneuro
{
  template <class S, class SMF>
  struct TransferMatrixUserTraits {
    using Solver = S;
    static const unsigned int dimension = S::Traits::dimension;
    using DenseRHSVector = typename Solver::Traits::RangeDOFVector;
    using SparseRHSVector = SparseVectorContainer<typename DenseRHSVector::ContainerIndex,
                                                  typename DenseRHSVector::ElementType>;
    using CoordinateFieldType = typename S::Traits::CoordinateFieldType;
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
    using DomainField = typename Solver::Traits::DomainDOFVector::field_type;
  };

  template <class S, class SMF>
  class TransferMatrixUser
  {
  public:
    using Traits = TransferMatrixUserTraits<S, SMF>;

    explicit TransferMatrixUser(std::shared_ptr<const typename Traits::Solver> solver)
        : solver_(solver)
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
            *solver_, config, solverConfig);
      } else {
        denseSourceModel_ = SMF::template createDense<typename Traits::DenseRHSVector>(
            *solver_, config, solverConfig);
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

    void postProcessPotential(const std::vector<ProjectedElectrode<typename S::Traits::GridView>>& projectedElectrodes,
                              std::vector<typename Traits::DomainField>& potential)
    {
      if (density_ == VectorDensity::sparse) {
        sparseSourceModel_->postProcessSolution(projectedElectrodes, potential);
      } else {
        denseSourceModel_->postProcessSolution(projectedElectrodes, potential);
      }
    }

    void postProcessMEG(const std::vector<typename Traits::Coordinate>& coils,
                        const std::vector<std::vector<typename Traits::Coordinate>>& projections,
                        std::vector<typename Traits::DomainField>& fluxes)
    {
      if(density_ == VectorDensity::sparse) {
        sparseSourceModel_->postProcessMEG(coils, projections, fluxes);
      }
      else {
        denseSourceModel_->postProcessMEG(coils, projections, fluxes);
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
      sparseSourceModel_->assembleRightHandSide(rhs);

      const auto blockSize = Traits::DenseRHSVector::block_type::dimension;

      std::vector<typename Traits::DomainField> output;
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
      denseSourceModel_->assembleRightHandSide(*denseRHSVector_);

      return matrix_dense_vector_product(transferMatrix,
                                         Dune::PDELab::Backend::native(*denseRHSVector_));
    }

  private:
    std::shared_ptr<const typename Traits::Solver> solver_;
    VectorDensity density_;
    std::shared_ptr<SourceModelInterface<typename S::Traits::GridView, typename Traits::DomainField, Traits::dimension,
                                         typename Traits::SparseRHSVector>>
        sparseSourceModel_;
    std::shared_ptr<SourceModelInterface<typename S::Traits::GridView, typename Traits::DomainField, Traits::dimension,
                                         typename Traits::DenseRHSVector>>
        denseSourceModel_;
    mutable std::shared_ptr<typename Traits::DenseRHSVector> denseRHSVector_;
  };
}

#endif // DUNEURO_TRANSFER_MATRIX_USER_HH
