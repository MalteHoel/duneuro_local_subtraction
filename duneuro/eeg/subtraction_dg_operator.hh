#ifndef DUNEURO_SUBTRACTION_DG_OPERATOR_HH
#define DUNEURO_SUBTRACTION_DG_OPERATOR_HH

/*
 * subtraction_dg_operator.hh
 *
 *	The local operator class which will be used by the grid operator to assemble the matrix
 *	and the residual.
 *	The class uses parts of the ConvectionDiffusionDG Operator implementation of PDELab.
 *	Some methods of the PDELab implementation had to be overwritten in order to adapt to our
 *	problem's formulation. For the most part this applies to the lambda_-methods which only
 *	depend on the test functions. They are implemented in the SubtractionDGLambda class.
 *	Additionally to the already implemented functions in the ConvectionDiffusionDG Operator
 *	we need a lambda_skeleton method for the boundary terms arising in the DG formulation of
 *	our problem.
 *
 *  Created on: Apr 25, 2013
 *      Author: jakob
 */

#include <duneuro/common/convection_diffusion_dg_operator.hh>
#include <duneuro/eeg/subtraction_dg_lambda.hh>

namespace duneuro
{
  /**** class definition of the operator ****/
  template <typename PROBLEMDATA, typename EdgeNormProvider>
  class SubtractionDG : public ConvectionDiffusion_DG_LocalOperator<PROBLEMDATA, EdgeNormProvider>,
                        public SubtractionDGLambda<PROBLEMDATA>
  {
  public:
    /*** flags that tell the grid operator what to do ***/
    enum { doAlphaBoundary = false };
    enum { doAlphaVolume = true };
    enum { doLambdaBoundary = true };
    enum { doLambdaVolume = true };
    enum { doAlphaSkeleton = true };
    enum { doLambdaSkeleton = true };

    using SubtractionDGLambda<PROBLEMDATA>::lambda_volume;
    // using SubtractionDGLambda<PROBLEMDATA>::lambda_boundary;
    using SubtractionDGLambda<PROBLEMDATA>::lambda_skeleton;

    /*** Constructor ***/
    SubtractionDG(
        PROBLEMDATA& problem_, const EdgeNormProvider& edgenormprovider_,
        const ConvectionDiffusion_DG_Scheme::Type scheme_ = ConvectionDiffusion_DG_Scheme::SIPG,
        const ConvectionDiffusion_DG_Weights::Type weights_ =
            ConvectionDiffusion_DG_Weights::weightsOn,
        double alpha_ = 0.0, unsigned int intorderadd_ = 0, unsigned int intorderadd_lb_ = 0)
        : ConvectionDiffusion_DG_LocalOperator<PROBLEMDATA, EdgeNormProvider>(
              problem_, edgenormprovider_, scheme_, weights_, alpha_, intorderadd_)
        , SubtractionDGLambda<PROBLEMDATA>(problem_, intorderadd_, intorderadd_lb_)
    {
    }
  };
}

#endif // DUNEURO_SUBTRACTION_DG_OPERATOR_HH
