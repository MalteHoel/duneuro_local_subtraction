#ifndef DUNEURO_CONFORMING_TRANSFER_MATRIX_USER_HH
#define DUNEURO_CONFORMING_TRANSFER_MATRIX_USER_HH

#include <duneuro/common/dg_solver.hh>
#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/sparse_vector_container.hh>
#include <duneuro/common/transfer_matrix.hh>
#include <duneuro/common/vector_density.hh>
#include <duneuro/eeg/dg_source_model_factory.hh>

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
    using Dipole = Dipole<CoordinateField, dimension>;
    using DomainField = typename EEGForwardSolver::Traits::DomainDOFVector::field_type;
    using TransferMatrix = ISTLTransferMatrix<DomainField>;
  };

  template <class S, class SMF>
  class ConformingTransferMatrixUser
  {
  public:
    using Traits = ConformingTransferMatrixUserTraits<S>;

    ConformingTransferMatrixUser(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                                 std::shared_ptr<typename Traits::TransferMatrix> transferMatrix,
                                 const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor)
        , transferMatrix_(transferMatrix)
        , solver_(volumeConductor_, config)
        , density_(source_model_default_density(config.sub("source_model")))
        , sourceModelDense_(SMF::template createDense<typename Traits::DenseRHSVector>(
              volumeConductor, solver_, config.sub("source_model")))
        , sourceModelSparse_(SMF::template createSparse<typename Traits::SparseRHSVector>(
              volumeConductor, solver_, config.sub("source_model")))
    {
    }

    std::vector<typename Traits::DomainField> solve(const typename Traits::Dipole& dipole) const
    {
      if (density_ == VectorDensity::sparse) {
        return solveSparse(dipole);
      } else {
        return solveDense(dipole);
      }
    }

    std::vector<typename Traits::DomainField>
    solveSparse(const typename Traits::Dipole& dipole) const
    {
      using SVC = typename Traits::SparseRHSVector;
      SVC rhs;
      sourceModelSparse_->assembleRightHandSide(dipole, rhs);

      const auto blockSize =
          Traits::EEGForwardSolver::Traits::FunctionSpace::GFS::Traits::Backend::blockSize;

      std::vector<typename Traits::DomainField> output(transferMatrix_->matrix().rows());
      if (blockSize == 1) {
        matrix_sparse_vector_product(transferMatrix_->matrix(), rhs, output,
                                     [](const typename SVC::Index& c) { return c[0]; });
      } else {
        matrix_sparse_vector_product(
            transferMatrix_->matrix(), rhs, output,
            [blockSize](const typename SVC::Index& c) { return c[1] * blockSize + c[0]; });
      }
      return output;
    }

    std::vector<typename Traits::DomainField>
    solveDense(const typename Traits::Dipole& dipole) const
    {
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(solver_, 0.0);
      }
      sourceModelDense_->assembleRightHandSide(dipole, *denseRHSVector_);

      std::vector<typename Traits::DomainField> output;
      output.reserve(transferMatrix_->matrix().rows());
      const auto blockSize =
          Traits::EEGForwardSolver::Traits::FunctionSpace::GFS::Traits::Backend::blockSize;
      for (std::size_t k = 0; k < transferMatrix_->matrix().rows(); ++k) {
        typename Traits::DomainField product = 0.0;
        for (std::size_t cb = 0; cb < denseRHSVector_->N(); ++cb) {
          for (std::size_t bi = 0; bi < blockSize; ++bi) {
            product +=
                transferMatrix_->matrix()[k][cb * blockSize + bi] * denseRHSVector_->block(cb)[bi];
          }
        }
        output.push_back(product);
      }
      return output;
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
    std::shared_ptr<typename Traits::TransferMatrix> transferMatrix_;
    typename Traits::EEGForwardSolver solver_;
    VectorDensity density_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateField, Traits::dimension,
                                         typename Traits::DenseRHSVector>>
        sourceModelDense_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateField, Traits::dimension,
                                         typename Traits::SparseRHSVector>>
        sourceModelSparse_;
    mutable std::shared_ptr<typename Traits::DenseRHSVector> denseRHSVector_;
  };
}

#endif // DUNEURO_CONFORMING_TRANSFER_MATRIX_USER_HH
