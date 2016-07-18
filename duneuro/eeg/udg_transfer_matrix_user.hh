#ifndef DUNEURO_UDG_TRANSFER_MATRIX_USER_HH
#define DUNEURO_UDG_TRANSFER_MATRIX_USER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <duneuro/common/dipole.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/make_dof_vector.hh>
#include <duneuro/common/matrix_utilities.hh>
#include <duneuro/common/sparse_vector_container.hh>
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
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGTransferMatrixUser
  {
  public:
    using Traits = UDGTransferMatrixUserTraits<ST, compartments, degree, DF, RF, JF>;

    UDGTransferMatrixUser(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                          std::shared_ptr<typename Traits::Solver> solver,
                          const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation)
        , solver_(solver)
        , density_(source_model_default_density(config.sub("source_model")))
        , sourceModelDense_(
              UDGSourceModelFactory::template createDense<compartments - 1,
                                                          typename Traits::DenseRHSVector>(
                  *solver_, config.sub("source_model")))
        , sourceModelSparse_(
              UDGSourceModelFactory::template createSparse<compartments - 1,
                                                           typename Traits::SparseRHSVector>(
                  *solver_, config.sub("source_model")))
    {
    }

    UDGTransferMatrixUser(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                          const Dune::ParameterTree& config)
        : UDGTransferMatrixUser(subTriangulation,
                                std::make_shared<typename Traits::Solver>(subTriangulation, config),
                                config)
    {
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
      sourceModelSparse_->assembleRightHandSide(dipole, rhs);

      const auto blockSize = Traits::Solver::Traits::FunctionSpace::blockSize;

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
    std::vector<typename Traits::DomainField>
    solveDense(const M& transferMatrix, const typename Traits::DipoleType& dipole) const
    {
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(*solver_, 0.0);
      }
      sourceModelDense_->assembleRightHandSide(dipole, *denseRHSVector_);

      return matrix_dense_vector_product(transferMatrix,
                                         Dune::PDELab::Backend::native(*denseRHSVector_));
    }

  private:
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::Solver> solver_;
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
