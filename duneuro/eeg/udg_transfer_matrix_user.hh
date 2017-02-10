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
    using Coordinate = Dune::FieldVector<CoordinateFieldType, dimension>;
    using DipoleType = Dipole<CoordinateFieldType, dimension>;
    using DomainField = typename Solver::Traits::DomainDOFVector::field_type;
    using ElementSearch = KDTreeElementSearch<typename ST::BaseT::GridView>;
  };

  template <class ST, int compartments, int degree, class DF = double, class RF = double,
            class JF = double>
  class UDGTransferMatrixUser
  {
  public:
    using Traits = UDGTransferMatrixUserTraits<ST, compartments, degree, DF, RF, JF>;

    UDGTransferMatrixUser(std::shared_ptr<typename Traits::SubTriangulation> subTriangulation,
                          std::shared_ptr<typename Traits::Solver> solver,
                          std::shared_ptr<typename Traits::ElementSearch> search,
                          const Dune::ParameterTree& config)
        : subTriangulation_(subTriangulation), solver_(solver), search_(search)
    {
    }

    void setSourceModel(const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      sparseSourceModel_.reset();
      denseSourceModel_.reset();
      density_ = source_model_default_density(config.sub("source_model"));
      if (density_ == VectorDensity::sparse) {
        sparseSourceModel_ =
            UDGSourceModelFactory::template createSparse<compartments - 1,
                                                         typename Traits::SparseRHSVector>(
                *solver_, subTriangulation_, search_, config);
      } else {
        denseSourceModel_ =
            UDGSourceModelFactory::template createDense<compartments - 1,
                                                        typename Traits::DenseRHSVector>(
                *solver_, subTriangulation_, search_, config);
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
      if (density_ == VectorDensity::sparse) {
        sparseSourceModel_->postProcessSolution(projectedElectrodes, potential);
      } else {
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
      sparseSourceModel_->assembleRightHandSide(rhs);

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
    std::shared_ptr<typename Traits::SubTriangulation> subTriangulation_;
    std::shared_ptr<typename Traits::Solver> solver_;
    std::shared_ptr<KDTreeElementSearch<typename ST::BaseT::GridView>> search_;
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

#endif // DUNEURO_UDG_TRANSFER_MATRIX_USER_HH
