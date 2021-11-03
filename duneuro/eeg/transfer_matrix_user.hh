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
#include <mutex>					// for std::mutex
#include <fstream>					// include for writing to files

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
   					 	      std::ofstream& rhs_assembly_ofstream,
                                                    std::ofstream& tm_multiplication_ofstream,
                                                    DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      std::vector<typename Traits::DomainField> result;
      if (density_ == VectorDensity::sparse) {
        dataTree.set("density", "sparse");
        result = solveSparse(transferMatrix, rhs_assembly_ofstream, tm_multiplication_ofstream);
      } else {
        dataTree.set("density", "dense");
        result = solveDense(transferMatrix, rhs_assembly_ofstream, tm_multiplication_ofstream);
      }
      dataTree.set("time", timer.elapsed());
      return result;
    }

    template <class M>
    std::vector<typename Traits::DomainField> solveSparse(const M& transferMatrix,
    							    std::ofstream& rhs_assembly_ofstream,
    							    std::ofstream& tm_multiplication_ofstream) const
    {
      using SVC = typename Traits::SparseRHSVector;
      SVC rhs;
      
      Dune::Timer timer(false);
      std::mutex write_to_file_mutex;
      
      timer.start();
      sparseSourceModel_->assembleRightHandSide(rhs);
      timer.stop();
      double time_rhs_assembly = timer.lastElapsed();

      {
        std::lock_guard<std::mutex> lock(write_to_file_mutex);
        rhs_assembly_ofstream << time_rhs_assembly << "\n";
      }


      const auto blockSize = Traits::DenseRHSVector::block_type::dimension;

      timer.start();
      std::vector<typename Traits::DomainField> output;

      auto product = (blockSize == 1) ? matrix_sparse_vector_product(transferMatrix, rhs,
                                            [](const typename SVC::Index& c) { return c[0]; }) 
                                   : matrix_sparse_vector_product(
            transferMatrix, rhs,
            [blockSize](const typename SVC::Index& c) { return c[1] * blockSize + c[0]; });

      timer.stop();
      double time_matrix_vector_product = timer.lastElapsed();
      
      {    
        std::lock_guard<std::mutex> lock(write_to_file_mutex);
        tm_multiplication_ofstream << time_matrix_vector_product << "\n";
      }
      
      return product;
    }

    template <class M>
    std::vector<typename Traits::DomainField> solveDense(const M& transferMatrix,
    							   std::ofstream& rhs_assembly_ofstream,
    							   std::ofstream& tm_multiplication_ofstream) const
    {
      Dune::Timer timer(false);
      std::mutex write_to_file_mutex;
    
      timer.start();
      if (!denseRHSVector_) {
        denseRHSVector_ = make_range_dof_vector(*solver_, 0.0);
      } else {
        *denseRHSVector_ = 0.0;
      }
      denseSourceModel_->assembleRightHandSide(*denseRHSVector_);
      timer.stop();
      double time_rhs_assembly = timer.lastElapsed();
      
      {
        std::lock_guard<std::mutex> lock(write_to_file_mutex);
        rhs_assembly_ofstream << time_rhs_assembly << "\n";
      }
      
      
      timer.start();
      auto product_vector = matrix_dense_vector_product(transferMatrix,
                                         Dune::PDELab::Backend::native(*denseRHSVector_));
      timer.stop();
      double time_matrix_vector_product = timer.lastElapsed();
      {
        std::lock_guard<std::mutex> lock(write_to_file_mutex);
        tm_multiplication_ofstream << time_matrix_vector_product << "\n";
      }
    
      return product_vector;
    }

  private:
    std::shared_ptr<const typename Traits::Solver> solver_;
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

#endif // DUNEURO_TRANSFER_MATRIX_USER_HH
