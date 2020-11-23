#ifndef DUNEURO_TDCS_EVALUATION_INTERFACE_HH
#define DUNEURO_TDCS_EVALUATION_INTERFACE_HH

#include <dune/udg/simpletpmctriangulation.hh>
#include <duneuro/udg/simpletpmc_domain.hh>

namespace duneuro
{
  template <typename GV, typename GFS>
  struct TDCSEvaluationInterface {
#if HAVE_DUNE_UDG
    using SubTriangulation = Dune::UDG::SimpleTpmcTriangulation<GV, GV>;
#endif
    using Element = typename GV::template Codim<0>::Entity;
    using LocalCoordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;
    using Coordinate = Dune::FieldVector<typename GV::ctype, GV::dimension>;

  public:
#if HAVE_DUNE_UDG
    virtual std::vector<std::vector<double>> evaluate(const std::vector<Coordinate>& positions, 
                                                      const DenseMatrix<double>& EvaluationMatrix, const SubTriangulation& subTriangulation) const = 0;
#endif
    virtual ~TDCSEvaluationInterface(){};
  };
}

#endif // DUNEURO_TDCS_EVALUATION_INTERFACE_HH