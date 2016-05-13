#ifndef DUNEURO_UDG_TRANSFER_MATRIX_USER_HH
#define DUNEURO_UDG_TRANSFER_MATRIX_USER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/sparse_vector_container.hh>
#include <duneuro/common/transfer_matrix.hh>
#include <duneuro/common/udg_solver.hh>
#include <duneuro/common/vector_density.hh>
#include <duneuro/eeg/udg_source_model_factory.hh>
#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class ST, int compartments, int degree, class DF, class RF, class JF>
  struct UDGTransferMatrixUserTraits {
    static const unsigned int dimension = ST::dim;
    using Solver = UDGSolver<ST, compartments, degree, DF, RF, JF>;
    using SubTriangulation = ST;
    using DenseRHSVector = typename Solver::Traits::RangeDOFVector;
    using SparseRHSVector = SparseVectorContainer<typename DenseRHSVector::ContainerIndex,
                                                  typename DenseRHSVector::ElementType>;
    using CoordinateFieldType = typename ST::ctype;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
    using DomainField = typename Solver::Traits::DomainDOFVector::field_type;
    using TransferMatrix = ISTLTransferMatrix<DomainField>;
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGTransferMatrixUser
  {
  public:
    using Traits = UDGTransferMatrixUserTraits<ST, compartments, degree, DF, RF, JF>;

    UDGTransferMatrixUser(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                          std::shared_ptr<typename Traits::TransferMatrix> transferMatrix,
                          const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , transferMatrix_(transferMatrix)
        , solver_(subTriangulation_, config)
        , density_(source_model_default_density(config.sub("source_model")))
        , sourceModelDense_(
              UDGSourceModelFactory::template createDense<compartments - 1,
                                                          typename Traits::DenseRHSVector>(
                  solver_, config.sub("source_model")))
        , sourceModelSparse_(
              UDGSourceModelFactory::template createSparse<compartments - 1,
                                                           typename Traits::SparseRHSVector>(
                  solver_, config.sub("source_model")))
    {
    }

    std::vector<typename Traits::DomainField> solve(const typename Traits::DipoleType& dipole,
                                                    DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      if (density_ == VectorDensity::sparse) {
        dataTree.set("density", "sparse");
        return solveSparse(dipole);
      } else {
        dataTree.set("density", "dense");
        return solveDense(dipole);
      }
      dataTree.set("time", timer.elapsed());
    }

    std::vector<typename Traits::DomainField>
    solveSparse(const typename Traits::DipoleType& dipole) const
    {
      using SVC = typename Traits::SparseRHSVector;
      SVC rhs;
      sourceModelSparse_->assembleRightHandSide(dipole, rhs);

      const auto blockSize = Traits::Solver::Traits::FunctionSpace::blockSize;

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
    solveDense(const typename Traits::DipoleType& dipole) const
    {
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(solver_, 0.0);
      }
      sourceModelDense_->assembleRightHandSide(dipole, *denseRHSVector_);

      std::vector<typename Traits::DomainField> output;
      output.reserve(transferMatrix_->matrix().rows());
      const auto blockSize = Traits::Solver::Traits::FunctionSpace::blockSize;
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
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::TransferMatrix> transferMatrix_;
    typename Traits::Solver solver_;
    VectorDensity density_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::DenseRHSVector>>
        sourceModelDense_;
    std::shared_ptr<SourceModelInterface<typename Traits::CoordinateFieldType, Traits::dimension,
                                         typename Traits::SparseRHSVector>>
        sourceModelSparse_;
    mutable std::shared_ptr<typename Traits::DenseRHSVector> denseRHSVector_;
  };
}

#endif // DUNEURO_UDG_TRANSFER_MATRIX_USER_HH
