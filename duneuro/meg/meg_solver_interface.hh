#ifndef DUNEURO_MEG_SOLVER_INTERFACE_HH
#define DUNEURO_MEG_SOLVER_INTERFACE_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/backend/interface.hh>

#include <duneuro/io/vtk_writer.hh>

namespace duneuro
{
  template <class VC, class V>
  class MEGSolverInterface
  {
  public:
    using DomainType = Dune::FieldVector<typename VC::ctype, VC::dim>;

    virtual void bind(const std::vector<DomainType>& sensor,
                      const std::vector<std::vector<DomainType>>& projection) = 0;
    virtual void bind(const V& eegSolution) = 0;
    virtual typename Dune::FieldTraits<Dune::PDELab::Backend::Native<V>>::field_type
    solve(std::size_t coil, std::size_t projection) const = 0;
    virtual void assembleTransferMatrixRHS(std::size_t coil, std::size_t projection,
                                           V& rhs) const = 0;
    virtual void addFluxToVTKWriter(VTKWriter<VC>& writer) const = 0;
    virtual std::size_t numberOfCoils() const = 0;
    virtual std::size_t numberOfProjections(std::size_t coil) const = 0;
    virtual ~MEGSolverInterface()
    {
    }
  };
}

#endif // DUNEURO_MEG_SOLVER_INTERFACE_HH
