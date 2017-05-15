#ifndef DUNEURO_SUBTRACTION_DG_OPERATOR_HH
#define DUNEURO_SUBTRACTION_DG_OPERATOR_HH

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/eeg/subtraction_dg_lambda.hh>

namespace duneuro
{
  enum class SubtractionContinuityType { continuous, discontinuous };

  /**** class definition of the operator ****/
  template <typename PROBLEMDATA, typename EdgeNormProvider,
            SubtractionContinuityType continuityType = SubtractionContinuityType::discontinuous>
  class SubtractionDG : public SubtractionDGLambda<PROBLEMDATA>,
                        public Dune::PDELab::FullVolumePattern,
                        public Dune::PDELab::FullSkeletonPattern,
                        public Dune::PDELab::LocalOperatorDefaultFlags
  {
  public:
    /*** flags that tell the grid operator what to do ***/
    enum { doLambdaBoundary = true };
    enum { doLambdaVolume = true };
    enum { doLambdaSkeleton = continuityType == SubtractionContinuityType::discontinuous };

    using SubtractionDGLambda<PROBLEMDATA>::lambda_volume;
    using SubtractionDGLambda<PROBLEMDATA>::lambda_boundary;
    using SubtractionDGLambda<PROBLEMDATA>::lambda_skeleton;

    /*** Constructor ***/
    SubtractionDG(PROBLEMDATA& problem_, const ConvectionDiffusion_DG_Weights::Type weights_ =
                                             ConvectionDiffusion_DG_Weights::weightsOn,
                  unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0)
        : SubtractionDGLambda<PROBLEMDATA>(problem_, intorderadd_, intorderadd_lb_, weights_)
    {
    }
  };
}

#endif // DUNEURO_SUBTRACTION_DG_OPERATOR_HH
