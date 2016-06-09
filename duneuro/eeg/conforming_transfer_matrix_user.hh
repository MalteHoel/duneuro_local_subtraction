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
    using DipoleType = Dipole<CoordinateField, dimension>;
    using DomainField = typename EEGForwardSolver::Traits::DomainDOFVector::field_type;
    using TransferMatrix = ISTLTransferMatrix<DomainField>;
  };

  template <class S, class SMF>
  class ConformingTransferMatrixUser
  {
  public:
    using Traits = ConformingTransferMatrixUserTraits<S>;

    ConformingTransferMatrixUser(std::shared_ptr<typename Traits::VolumeConductor> volumeConductor,
                                 const Dune::ParameterTree& config)
        : volumeConductor_(volumeConductor)
        , solver_(volumeConductor_, config)
        , density_(source_model_default_density(config.sub("source_model")))
    {
      if (density_ == VectorDensity::dense) {
        sourceModelDense_ = SMF::template createDense<typename Traits::DenseRHSVector>(
            volumeConductor, solver_, config.sub("source_model"));
      } else {
        sourceModelSparse_ = SMF::template createSparse<typename Traits::SparseRHSVector>(
            volumeConductor, solver_, config.sub("source_model"));
      }
    }

    template <class M>
    std::vector<typename Traits::DomainField> solve(const M& transferMatrix,
                                                    const typename Traits::DipoleType& dipole,
                                                    DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      if (density_ == VectorDensity::sparse) {
        dataTree.set("density", "sparse");
        return solveSparse(transferMatrix, dipole);
      } else {
        dataTree.set("density", "dense");
        return solveDense(transferMatrix, dipole);
      }
      dataTree.set("time", timer.elapsed());
    }

    template <class M>
    std::vector<typename Traits::DomainField>
    solveSparse(const M& transferMatrix, const typename Traits::DipoleType& dipole) const
    {
      using SVC = typename Traits::SparseRHSVector;
      SVC rhs;
      assert(sourceModelSparse_);
      sourceModelSparse_->assembleRightHandSide(dipole, rhs);

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
    std::vector<typename Traits::DomainField>
    solveDense(const M& transferMatrix, const typename Traits::DipoleType& dipole) const
    {
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(solver_, 0.0);
      }
      assert(sourceModelDense_);
      sourceModelDense_->assembleRightHandSide(dipole, *denseRHSVector_);

      return matrix_dense_vector_product(transferMatrix,
                                         Dune::PDELab::Backend::native(*denseRHSVector_));
    }

  private:
    std::shared_ptr<typename Traits::VolumeConductor> volumeConductor_;
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
