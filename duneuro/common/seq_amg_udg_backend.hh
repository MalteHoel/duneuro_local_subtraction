// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_SEQ_AMG_UDG_BACKEND_HH
#define DUNEURO_SEQ_AMG_UDG_BACKEND_HH

#include <dune/common/parametertree.hh>
#include <dune/common/power.hh>

#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/overlappingschwarz.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/ovlpistlsolverbackend.hh>
#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <duneuro/common/cg_to_udg_prolongation.hh>

namespace duneuro
{
  namespace AMG4UDGDetail
  {
    template <class GFS, class ST, class E>
    std::set<std::size_t> extractBlocks(const GFS& gfs, const ST& st, const E& element)
    {
      using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
      LFS lfs(gfs);
      Dune::PDELab::LFSIndexCache<LFS> cache(lfs);
      lfs.bind(element);
      cache.update();
      std::set<std::size_t> result;
      for (unsigned int i = 0; i < cache.size(); ++i) {
        result.insert(cache.containerIndex(i).back());
      }
      return result;
    }

    template <class GFS, class ST>
    std::vector<std::set<std::size_t>> extractCouplings(const GFS& gfs, const ST& st)
    {
      Dune::SingleCodimSingleGeomTypeMapper<typename GFS::Traits::GridView, 0> elementMapper(
          gfs.gridView());
      std::vector<std::set<std::size_t>> result;
      std::vector<bool> considered(elementMapper.size());
      for (const auto& element : Dune::elements(gfs.gridView())) {
        auto insideIndex = elementMapper.index(element);
        auto current = extractBlocks(gfs, st, element);
        // if we have more than one subdomain within element
        if (current.size() > 1) {
          // find all intersections with neighbor containing more than one subdomain
          // note that this is not exact
          for (const auto& it : Dune::intersections(gfs.gridView(), element)) {
            if (it.neighbor()) {
              auto outsideIndex = elementMapper.index(it.outside());
              if (insideIndex < outsideIndex) {
                auto out = extractBlocks(gfs, st, it.outside());
                if (out.size() > 1) {
                  out.insert(current.begin(), current.end());
                  result.push_back(out);
                }
              }
            }
          }
          // if its exactly one subdomain
        } else if (current.size() == 1) {
          result.push_back(current);
        }
      }
      return result;
    }

    template <class T, int N, int M, class P>
    void assertAll(const Dune::BCRSMatrix<Dune::FieldMatrix<T, N, M>>& m, const std::string& name,
                   P predicate)
    {
      for (auto row = m.begin(); row != m.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
          for (int r = 0; r < N; ++r) {
            for (int c = 0; c < M; ++c) {
              if (!predicate((*col)[r][c])) {
                std::cout << "predicate \"" << name << "\" not fullfilled at block (" << row.index()
                          << "," << col.index() << ") entry (" << r << "," << c
                          << "): " << (*col)[r][c] << "\n";
              }
            }
          }
        }
      }
    }

    template <class T, int N, int M>
    void write(const std::string& filename,
               const Dune::BCRSMatrix<Dune::FieldMatrix<T, N, M>>& matrix)
    {
      std::ofstream out(filename);
      Dune::printmatrix(out, matrix, "jacobian", "r");
    }
  }
  /** Sequential solver backend for using AMG for DG in PDELab

      The template parameters are:
      DGGO       GridOperator for DG discretization, allows access to matrix, vector and grid
     function space
      CGGFS      grid function space for CG subspace
      DGPrec     preconditioner for DG problem
      Solver     solver to be used on the complete problem

  */
  template <class DGGO, class UDGGFS, class ST, template <class, class, class, int> class DGPrec,
            template <class> class Solver>
  class ISTLBackend_SEQ_AMG_4_UDG : public Dune::PDELab::LinearResultStorage
  {
    // DG grid function space
    typedef typename DGGO::Traits::TrialGridFunctionSpace GFS;

    // vectors and matrices on DG level
    typedef typename DGGO::Traits::Jacobian M; // wrapped istl DG matrix
    typedef typename DGGO::Traits::Domain V; // wrapped istl DG vector
    typedef Dune::PDELab::Backend::Native<M> Matrix; // istl DG matrix
    typedef Dune::PDELab::Backend::Native<V> Vector; // istl DG vector
    typedef typename Vector::field_type field_type;

    // prolongation matrix
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<field_type, M::block_type::rows, 1>> P;
    typedef Dune::BlockVector<Dune::FieldVector<field_type, 1>> CGVector;

    // CG subspace matrix
    typedef typename Dune::TransposedMatMultMatResult<P, Matrix>::type PTADG;
    typedef typename Dune::MatMultMatResult<PTADG, P>::type ACG; // istl coarse space matrix
    typedef ACG CGMatrix; // another name

    // AMG in CG-subspace
    typedef Dune::MatrixAdapter<CGMatrix, CGVector, CGVector> CGOperator;
    typedef Dune::SeqSSOR<CGMatrix, CGVector, CGVector, 1> Smoother;
    // typedef Dune::SeqOverlappingSchwarz<CGMatrix, CGVector,
    // Dune::SymmetricMultiplicativeSchwarzMode> Smoother;
    typedef Dune::Amg::AMG<CGOperator, CGVector, Smoother> AMG;
    typedef Dune::Amg::Parameters Parameters;

    const UDGGFS& udggfs;
    const ST& st;

    int verbose;
    unsigned int cgSmootherIterations;
    field_type cgSmootherRelaxation;
    unsigned int gamma;
    unsigned int dgSmootherIterations;
    field_type dgSmootherRelaxation;
    unsigned int solverIterations;
    bool additive;
    std::string smoother;

    bool reuse;
    bool firstapply;

    std::shared_ptr<CGOperator> cgop;
    std::shared_ptr<AMG> amg;
    std::shared_ptr<P> pmatrix;
    ACG acg; // CG-subspace matrix

  public:
    ISTLBackend_SEQ_AMG_4_UDG(const UDGGFS& udggfs_, const ST& st_,
                              const Dune::ParameterTree& params)
        : udggfs(udggfs_)
        , st(st_)
        , verbose(params.get<int>("verbose", 0))
        , cgSmootherIterations(params.get<unsigned int>("cg_smoother_iterations", 1))
        , cgSmootherRelaxation(params.get<field_type>("cg_smoother_relaxation", 1.0))
        , gamma(params.get<unsigned int>("gamma", 1))
        , dgSmootherIterations(params.get<unsigned int>("dg_smoother_iterations", 1))
        , dgSmootherRelaxation(params.get<field_type>("dg_smoother_relaxation", 1.0))
        , solverIterations(params.get<unsigned int>("max_iterations", 5000))
        , additive(params.get<bool>("additive", false))
        , smoother(params.get<std::string>("smoother"))
        , reuse(true)
        , firstapply(true)
    {
    }

    void setParameters(const Dune::ParameterTree& config)
    {
      if (config.hasKey("cg_smoother_iterations"))
        cgSmootherIterations = config.get<unsigned int>("cg_smoother_iterations");
      if (config.hasKey("cg_smoother_relaxation"))
        cgSmootherRelaxation = config.get<field_type>("cg_smoother_relaxation");
      if (config.hasKey("gamma"))
        gamma = config.get<unsigned int>("gamma");
      if (config.hasKey("dg_smoother_iterations"))
        dgSmootherIterations = config.get<unsigned int>("dg_smoother_iterations");
      if (config.hasKey("dg_smoother_relaxation"))
        dgSmootherRelaxation = config.get<field_type>("dg_smoother_relaxation");
      if (config.hasKey("max_iterations"))
        solverIterations = config.get<unsigned int>("max_iterations");
      if (config.hasKey("additive"))
        additive = config.get<bool>("additive");
      if (config.hasKey("verbose"))
        verbose = config.get<int>("verbose");
      if (config.hasKey("smoother"))
        smoother = config.get<std::string>("smoother");
    }

    /*! \brief compute global norm of a vector

      \param[in] v the given vector
    */
    typename V::ElementType norm(const V& v) const
    {
      return Dune::PDELab::Backend::native(v).two_norm();
    }

    //! Set whether the AMG should be reused again during call to apply().
    void setReuse(bool reuse_)
    {
      reuse = reuse_;
    }

    //! Return whether the AMG is reused during call to apply()
    bool getReuse() const
    {
      return reuse;
    }

    /*! \brief solve the given linear system

      \param[in] A the given matrix
      \param[out] z the solution vector to be computed
      \param[in] r right hand side
      \param[in] reduction to be achieved
    */
    void apply(M& A, V& z, V& r, typename V::ElementType reduction)
    {
      using Dune::PDELab::Backend::native;
      // do triple matrix product ACG = P^T ADG P
      Dune::Timer watch;
      watch.reset();
      // only do triple matrix product if the matrix changes
      double triple_product_time = 0.0;
      // no need to set acg here back to zero, this is done in matMultmat
      if (reuse == false || firstapply == true) {
        // pmatrix = constructCG2UDGProlongation(native(A));
        pmatrix = compute_cg_to_udg_prolongation(udggfs, st);
        //auto pmatrix_new = compute_cg_to_udg_prolongation_new(udggfs, st);
        //*pmatrix_new -= *pmatrix;
        //std::cout << "difference: " << pmatrix_new->frobenius_norm() << "\n";
        // Dune::writeMatrixToMatlab(*pmatrix, "pmatrix.mat");
        // Dune::writeMatrixToMatlab(native(A), "A.mat");
        std::cout << "A size: " << native(A).N() << "x" << native(A).M() << std::endl;
        std::cout << "pmatrix size: " << pmatrix->N() << "x" << pmatrix->M() << std::endl;
        PTADG ptadg;
        Dune::transposeMatMultMat(ptadg, *pmatrix, native(A));
        std::cout << "ptadg size: " << ptadg.N() << "x" << ptadg.M() << std::endl;
        Dune::matMultMat(acg, ptadg, *pmatrix);
        std::cout << "acg size: " << acg.N() << "x" << acg.M() << std::endl;
        // AMG4UDGDetail::write("acg.txt", acg);
        triple_product_time = watch.elapsed();
        if (verbose > 0)
          std::cout << "=== triple matrix product " << triple_product_time << " s" << std::endl;
        std::cout << "first block of acg:\n" << acg[0][0] << std::endl;
        // std::ofstream fstr("tmpmat.txt");
        // Dune::printmatrix(fstr, acg, "triple product matrix", "row", 10, 2);
        AMG4UDGDetail::assertAll(acg, "not nan", [](double v) { return !std::isnan(v); });
        // Dune::writeMatrixToMatlab(acg, "acg.mat");
      } else if (verbose > 0)
        std::cout << "=== reuse CG matrix, SKIPPING triple matrix product " << std::endl;

      // set up AMG solver for the CG subspace
      typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
      SmootherArgs smootherArgs;
      smootherArgs.iterations = cgSmootherIterations;
      smootherArgs.relaxationFactor = cgSmootherRelaxation;
      typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<CGMatrix,
                                                                        Dune::Amg::FirstDiagonal>>
          Criterion;
      Parameters amg_parameters(15, 2000);
      Criterion criterion(amg_parameters);
      amg_parameters.setGamma(gamma);
      amg_parameters.setAdditive(additive);
      watch.reset();

      // only construct a new AMG for the CG-subspace if the matrix changes
      double amg_setup_time = 0.0;
      if (reuse == false || firstapply == true) {
        cgop = std::make_shared<CGOperator>(acg);
        amg = std::make_shared<AMG>(*cgop, criterion, smootherArgs);
        firstapply = false;
        amg_setup_time = watch.elapsed();
        if (verbose > 0)
          std::cout << "=== AMG setup " << amg_setup_time << " s" << std::endl;
      } else if (verbose > 0)
        std::cout << "=== reuse CG matrix, SKIPPING AMG setup " << std::endl;

      // set up hybrid DG/CG preconditioner
      Dune::MatrixAdapter<Matrix, Vector, Vector> op(native(A));
      Dune::InverseOperatorResult stat;
      if (smoother == "schwarz") {
        using DGPreconditioner =
            Dune::SeqOverlappingSchwarz<Matrix, Vector, Dune::SymmetricMultiplicativeSchwarzMode>;
        auto domains = AMG4UDGDetail::extractCouplings(udggfs, st);
        DGPreconditioner dgprec(native(A), domains, dgSmootherRelaxation);

        typedef Dune::PDELab::SeqDGAMGPrec<Matrix, DGPreconditioner, AMG, P> HybridPrec;
        HybridPrec hybridprec(native(A), dgprec, *amg, *pmatrix, dgSmootherIterations,
                              dgSmootherIterations);

        // set up solver
        Solver<Vector> solver(op, hybridprec, reduction, solverIterations, verbose);

        // solve
        watch.reset();
        solver.apply(native(z), native(r), stat);
      } else if (smoother == "default") {
        DGPrec<Matrix, Vector, Vector, 1> dgprec(native(A), 1, 1);
        typedef Dune::PDELab::SeqDGAMGPrec<Matrix, DGPrec<Matrix, Vector, Vector, 1>, AMG, P>
            HybridPrec;
        HybridPrec hybridprec(native(A), dgprec, *amg, *pmatrix, dgSmootherIterations,
                              dgSmootherIterations);

        // set up solver
        Solver<Vector> solver(op, hybridprec, reduction, solverIterations, verbose);

        // solve
        watch.reset();
        solver.apply(native(z), native(r), stat);
      } else {
        DUNE_THROW(Dune::Exception, "unknown smoother type \"" << smoother << "\"");
      }
      double amg_solve_time = watch.elapsed();
      if (verbose > 0)
        std::cout << "=== Hybrid total solve time "
                  << amg_solve_time + amg_setup_time + triple_product_time << " s" << std::endl;
      res.converged = stat.converged;
      res.iterations = stat.iterations;
      res.elapsed = amg_solve_time + amg_setup_time + triple_product_time;
      res.reduction = stat.reduction;
      res.conv_rate = stat.conv_rate;
    }
  };
}
#endif
