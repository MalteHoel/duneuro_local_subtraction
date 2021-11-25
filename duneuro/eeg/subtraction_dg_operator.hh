#ifndef DUNEURO_SUBTRACTION_DG_OPERATOR_HH
#define DUNEURO_SUBTRACTION_DG_OPERATOR_HH

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/eeg/subtraction_dg_lambda.hh>
#include <duneuro/common/flags.hh>

namespace duneuro
{

  /**** class definition of the operator ****/
  template <typename PROBLEMDATA, typename EdgeNormProvider, typename PenaltyFluxWeighting,
            ContinuityType continuityType = ContinuityType::discontinuous>
  class SubtractionDG : public SubtractionDGLambda<PROBLEMDATA, PenaltyFluxWeighting>,
                        public Dune::PDELab::FullVolumePattern,
                        public Dune::PDELab::FullSkeletonPattern,
                        public Dune::PDELab::LocalOperatorDefaultFlags
  {
  public:
    /*** flags that tell the grid operator what to do ***/
    enum { doLambdaBoundary = true };
    enum { doLambdaVolume = true };
    enum { doLambdaSkeleton = continuityType == ContinuityType::discontinuous };

    using SubtractionDGLambda<PROBLEMDATA, PenaltyFluxWeighting>::lambda_volume;
    using SubtractionDGLambda<PROBLEMDATA, PenaltyFluxWeighting>::lambda_boundary;
    using SubtractionDGLambda<PROBLEMDATA, PenaltyFluxWeighting>::lambda_skeleton;

    /*** Constructor ***/
    SubtractionDG(PROBLEMDATA& problem_, const PenaltyFluxWeighting& weighting_,
                  unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0)
        : SubtractionDGLambda<PROBLEMDATA, PenaltyFluxWeighting>(problem_, weighting_, intorderadd_,
                                                                 intorderadd_lb_)
    {
    }
  };
}

#endif // DUNEURO_SUBTRACTION_DG_OPERATOR_HH
