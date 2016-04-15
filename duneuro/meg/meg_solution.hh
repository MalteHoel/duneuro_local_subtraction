#ifndef DUNEURO_MEGSOLUTION_HH
#define DUNEURO_MEGSOLUTION_HH

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>

#include <duneuro/io/data_tree.hh>
#include <duneuro/meg/meg_local_operator.hh>
#include <duneuro/meg/meg_transfer_matrix_rhs.hh>

#include <vector>

namespace duneuro
{
  template <class VC, class FS, class DF>
  class MEGSolution
  {
  public:
    using DOF = typename FS::DOF;
    using Coordinate = Dune::FieldVector<typename VC::ctype, VC::dim>;
    using LOP = MEGLocalOperator<VC>;
    using Assembler =
        Dune::PDELab::GalerkinGlobalAssembler<FS, LOP, Dune::SolverCategory::sequential>;

    MEGSolution(std::shared_ptr<VC> volumeConductor, const FS& fs,
                const std::vector<Coordinate>& coils,
                const std::vector<std::vector<Coordinate>>& projections,
                const Dune::ParameterTree& config, DataTree dataTree = DataTree())
        : config_(config), smatrix_(coils.size())
    {
      Dune::Timer timer;
      DOF smatrixrow(fs.getGFS(), 0.0);

      MEGTransferMatrixRHS<VC, FS> rhsAssembler(volumeConductor, &fs, config);
      for (unsigned int i = 0; i < coils.size(); ++i) {
        std::cout << "\rworking on coil " << i << "/" << coils.size();
        std::flush(std::cout);
        smatrix_[i].reserve(projections[i].size());
        for (const auto& p : projections[i]) {
          smatrixrow = 0.0;
          rhsAssembler.assembleRightHandSide(coils[i], p, smatrixrow);
          smatrix_[i].push_back(smatrixrow);
        }
      }
      std::cout << "\n";
      timer.stop();
      dataTree.set("time", timer.elapsed());
    }

    template <class I>
    void evaluate(const DOF& x, I output, DataTree dataTree = DataTree()) const
    {
      Dune::Timer timer;
      for (unsigned int coilIndex = 0; coilIndex < smatrix_.size(); ++coilIndex) {
        std::vector<typename DOF::field_type> result;
        result.reserve(smatrix_[coilIndex].size());
        for (unsigned int i = 0; i < smatrix_[coilIndex].size(); ++i) {
          result.push_back(smatrix_[coilIndex][i] * x);
        }
        *output++ = result;
      }
      dataTree.set("time", timer.elapsed());
    }

    std::vector<std::vector<typename DOF::field_type>>
    evaluate(const DOF& x, DataTree dataTree = DataTree()) const
    {
      std::vector<std::vector<typename DOF::field_type>> out;
      out.reserve(smatrix_.size());
      evaluate(x, std::back_inserter(out), dataTree);
      return out;
    }

  private:
    Dune::ParameterTree config_;
    std::vector<std::vector<DOF>> smatrix_;
  };
}

#endif // DUNEURO_MEGSOLUTION_HH
