#ifndef DUNEURO_STATIONARYLINEARPROBLEM_HH
#define DUNEURO_STATIONARYLINEARPROBLEM_HH

#include <iostream>
#include <memory>
#include <mutex>
#include <thread>

#include <dune/common/float_cmp.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/istl/bcrsmatrix.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/solver.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  struct UnsymmetricMatrixException : public Dune::Exception {
  };

  namespace TSSLPDetail
  {
    template <class T>
    struct SymmetryStatistics {
      T maximalAbsoluteDifference;
      std::pair<std::size_t, std::size_t> blockIndex;
      std::pair<std::size_t, std::size_t> localIndex;

      template <int N>
      explicit SymmetryStatistics(Dune::BCRSMatrix<Dune::FieldMatrix<T, N, N>>& m)
          : maximalAbsoluteDifference(0.0)
      {
        for (auto rit = m.begin(); rit != m.end(); ++rit) {
          for (auto cit = rit->begin(); cit != rit->end(); ++cit) {
            if (!m.exists(cit.index(), rit.index())) {
              DUNE_THROW(UnsymmetricMatrixException, "sparsity pattern is not symmetric. entry ("
                                                         << cit.index() << "," << rit.index()
                                                         << ") does not exist");
            } else {
              const auto& thism = *cit;
              const auto& otherm = m[cit.index()][rit.index()];
              // check that otherm == thism^t
              for (unsigned int r = 0; r < N; ++r) {
                for (unsigned int c = 0; c < N; ++c) {
                  auto diff = std::abs(thism[r][c] - otherm[c][r]);
                  if (diff > maximalAbsoluteDifference) {
                    maximalAbsoluteDifference = diff;
                    blockIndex = {rit.index(), cit.index()};
                    localIndex = {r, c};
                  }
                }
              }
            }
          }
        }
      }

      void print() const
      {
        std::cout << "SymmetryStatistics: maximal absolute difference: "
                  << maximalAbsoluteDifference << " at block (" << blockIndex.first << ","
                  << blockIndex.second << "), local index (" << localIndex.first << ","
                  << localIndex.second << ")" << std::endl;
      }
    };

    struct IllegalEntryException : public Dune::Exception {
    };

    template <class T, int N, class F>
    void assertEachEntry(const Dune::BCRSMatrix<Dune::FieldMatrix<T, N, N>>& m, F predicate)
    {
      for (auto rit = m.begin(); rit != m.end(); ++rit) {
        for (auto cit = rit->begin(); cit != rit->end(); ++cit) {
          for (int rb = 0; rb < N; ++rb) {
            for (int cb = 0; cb < N; ++cb) {
              if (!predicate((*cit)[rb][cb])) {
                DUNE_THROW(IllegalEntryException, "illegal entry found at block ("
                                                      << rit.index() << "," << cit.index()
                                                      << ") local index (" << rb << "," << cb
                                                      << ") : " << (*cit)[rb][cb]);
              }
            }
          }
        }
      }
    }

    template <class T, int N>
    void fixFirstDOF(Dune::BCRSMatrix<Dune::FieldMatrix<T, N, N>>& m, T value)
    {
      for (auto rit = m.begin(); rit != m.end(); ++rit) {
        if (rit.index() == 0) {
          for (auto cit = rit->begin(); cit != rit->end(); ++cit) {
            (*cit)[0][0] = (cit.index() == 0 ? value : 0.0);
            for (int c = 1; c < N; ++c) {
              (*cit)[0][c] = 0.0;
            }
            if (cit.index() == 0) {
              for (int r = 1; r < N; ++r) {
                m[rit.index()][0][r][0] = 0.0;
              }
            }
          }
        } else {
          // set first column to zero
          if (m.exists(rit.index(), 0)) {
            for (int r = 0; r < N; ++r) {
              m[rit.index()][0][r][0] = 0.0;
            }
          }
        }
      }
    }
  }

  //===============================================================
  // A class for solving linear stationary problems.
  // It assembles the matrix, computes the right hand side and
  // solves the problem.
  // This is only a first vanilla implementation which has to be improved.
  //===============================================================

  template <typename GO, typename DV, typename RV>
  class LinearProblemSolver
  {
    typedef typename GO::Traits::Jacobian M;
    typedef typename GO::Traits::TrialGridFunctionSpace TrialGridFunctionSpace;
    typedef GO GridOperator;

  public:
    typedef Dune::PDELab::StationaryLinearProblemSolverResult<double> Result;

    LinearProblemSolver(const GO& go, typename RV::ElementType reduction, bool fixDOF,
                        typename M::field_type fixedDOFEntry, int verbose = 1, bool debug = false)
        : _go(go)
        , _fixFirstDOF(fixDOF)
        , _fixedDOFEntry(fixedDOFEntry)
        , _verbose(verbose)
        , _debug(debug)
    {
    }

    //! Construct a StationaryLinearProblemSolver for the given objects and read parameters from
    //! a
    //! ParameterTree.
    /**
     * This constructor reads the parameter controlling its operation from a passed-in
     * ParameterTree
     * instead of requiring the user to specify all of them as individual constructor
     * parameters.
     * Currently the following parameters are read:
     *
     * Name                       | Default Value | Explanation
     * -------------------------- | ------------- | -----------
     * reduction                  |               | Required relative defect reduction
     * min_defect                 | 1e-99         | minimum absolute defect at which to stop
     * hanging_node_modifications | false         | perform required transformations for hanging
     * nodes
     * keep_matrix                | true          | keep matrix between calls to apply() (but
     * reassemble values every time)
     * verbosity                  | 1             | control amount of debug output
     *
     * Apart from reduction, all parameters have a default value and are optional.
     * The actual reduction for a call to apply() is calculated as r =
     * max(reduction,min_defect/start_defect),
     * where start defect is the norm of the residual of x.
     */
    LinearProblemSolver(const GO& go, const Dune::ParameterTree& params)
        : _go(go)
        , _fixFirstDOF(params.get<bool>("fixDOF", true))
        , _fixedDOFEntry(params.get<typename M::field_type>("fixedDOFEntry", 1.0))
        , _verbose(params.get<int>("verbosity", 1))
        , _debug(params.get<bool>("debug", false))
    {
    }

    template <class LS>
    void apply(LS& ls, DV& x, const Dune::ParameterTree& config, DataTree dataTree = DataTree())
    {
      DV tmp(_go.trialGridFunctionSpace(), 0.0);
      RV rhs(_go.testGridFunctionSpace(), 0.0);
      _go.residual(tmp, rhs);
      rhs *= -1.0;
      apply(ls, x, rhs, config, dataTree);
    }

    template <class LS>
    void apply(LS& ls, DV& x, const RV& rightHandSide, const Dune::ParameterTree& config,
               DataTree dataTree = DataTree())
    {
      Dune::Timer timer(false);
      {
        std::lock_guard<std::mutex> lock(_jacobian_mutex);
        if (!_jacobian) {
          timer.start();
          _jacobian = std::make_unique<M>(_go);
          timer.stop();
          if (_go.trialGridFunctionSpace().gridView().comm().rank() == 0 && _verbose >= 1)
            std::cout << "=== matrix setup (max) " << timer.lastElapsed() << " s" << std::endl;
          dataTree.set("time_matrix_setup", timer.lastElapsed());
          timer.start();
          (*_jacobian) = typename M::field_type(0.0);
          _go.jacobian(x, *_jacobian);
          if (_fixFirstDOF) {
            TSSLPDetail::fixFirstDOF(Dune::PDELab::Backend::native(*_jacobian), _fixedDOFEntry);
          }
          if (_debug) {
            TSSLPDetail::SymmetryStatistics<typename M::field_type> statistics(
                Dune::PDELab::Backend::native(*_jacobian));
            statistics.print();
            try {
              TSSLPDetail::assertEachEntry(Dune::PDELab::Backend::native(*_jacobian),
                                           [](typename M::field_type v) { return !std::isnan(v); });
            } catch (TSSLPDetail::IllegalEntryException& ex) {
              std::cout << "Illegal entry found:\n" << ex.what() << "\n";
            }
            std::cout << Dune::PDELab::Backend::native(*_jacobian)[0][0] << "\n";
          }
          timer.stop();
          dataTree.set("time_matrix_assembly", timer.lastElapsed());
        }
      }

      // transform rhs to discrete residuum
      RV r(rightHandSide);
      Dune::PDELab::Backend::native(*_jacobian)
          .mmv(Dune::PDELab::Backend::native(x), Dune::PDELab::Backend::native(r));
      r *= -1.0;
      // compute correction
      Dune::Timer solutionTimer;
      timer.start();
      DV z(_go.trialGridFunctionSpace(), 0.0);
      ls.apply(*_jacobian, z, r, config.get<typename RV::ElementType>(
                                     "reduction")); // solver makes right hand side consistent
      timer.stop();
      dataTree.set("iterations", ls.result().iterations);
      dataTree.set("reduction", ls.result().reduction);
      dataTree.set("conv_rate", ls.result().conv_rate);
      dataTree.set("time_solution", timer.lastElapsed());
      // and update
      timer.start();
      x -= z;
      timer.stop();
      dataTree.set("time", timer.elapsed());
    }

    //! Discard the stored Jacobian matrix.
    void discardMatrix()
    {
      std::lock_guard<std::mutex> lock(_jacobian_mutex);
      if (_jacobian)
        _jacobian.reset();
    }

  private:
    std::mutex _jacobian_mutex;
    const GO& _go;
    std::unique_ptr<M> _jacobian;
    bool _fixFirstDOF;
    typename M::field_type _fixedDOFEntry;
    Result _res;
    int _verbose;
    bool _debug;
  };
}

#endif
