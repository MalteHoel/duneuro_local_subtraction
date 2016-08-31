//$6    21.05.2003  Anwander A.  added empty setSystemMatrix()
//$5    11.03.2002  Matthias D.  Removed spelline error DoSetParamterValues();
//$4    14.01.2002    Matthias D.  Changed IP_SIM_SensorConfiguration_c to anEEGSensorConfiguration_c
//$3    03.08.2001    Matthias D.  Added DoSetParamterValues();
//$2    20.09.2000  Matthias D.  copied method from fweegforwardsimulator to IP_SIM_SimulatorEEGSpheres_c
//$1    28.08.2000  Matthias D.  created

#ifndef __IP_SIM_SIMULATOREEGSPHERES_C_H__
#define __IP_SIM_SIMULATOREEGSPHERES_C_H__

#include "CON_Vector_t.h"
#include "CON_Matrix_t.h"
#include "CON_Block_t.h"
#include "LIN_Auxil.h"

#include <math.h>

// IP_SIM_SimulatorEEGSpheres_c.h: interface for the IP_SIM_SimulatorEEGSpheres_c class.
//
//////////////////////////////////////////////////////////////////////



#include "IP_SIM_AbstractSimulatorEEGMEG_c.h"
#include "IP_SIM_SensorConfiguration_c.h"



class ANALYSIS_EXPORT IP_SIM_SimulatorEEGSpheres_c : public IP_SIM_AbstractSimulatorEEGMEG_c
{
public:
    IP_SIM_SimulatorEEGSpheres_c(anEEGSensorConfiguration_c& inSensorconfiguration,CON_Vector_t<double>& inRadii,CON_Vector_t<double>& inCenter,
                            CON_Vector_t<double>& inConductivities);
    virtual ~IP_SIM_SimulatorEEGSpheres_c();
    bool isValidSimulator();
    bool setRadii(CON_Vector_t<double>& inRadii);
    bool getRadii(CON_Vector_t<double>& outRadii);
    int  getNumberofSpheres();
    bool setCenter(CON_Vector_t<double>& inCenter);
    bool getCenter(CON_Vector_t<double>& outCenter);
    bool setConductivities(CON_Vector_t<double>& inConductivities);
    bool getConductivities(CON_Vector_t<double>& outConductivities);
    void computeGainMatrix(int NumberDip, const CON_Matrix_t<double>& inPos, const CON_Matrix_t<double>& inDir, CON_Matrix_t<double>& outSimData);
    void computeGainMatrix(const CON_Vector_t<double>& inParameters,CON_Matrix_t<double>& outSimData);

    void setSystemMatrix(const CON_Matrix_t<double>& SystemMatrix) {};
    void setSystemMatrix(const CON_Matrix_t<float>& SystemMatrix) {};

protected:
    void   computeEEGFourLayerSphere(double cta, double fi, double r0,double epslon[4], double eta[4], double r[4], CON_Vector_t<double>& pdip);
    double plgndr(int l, int m, double x);   // Compute Result of Legendre Series

    bool DoSetParameterValues(); // Transfers parameter values to variables used for calculations

protected:
anEEGSensorConfiguration_c& m_SensorConfiguration;
CON_Vector_t<double>&            m_Radii;
CON_Vector_t<double>&            m_Center;
CON_Vector_t<double>&            m_Conductivities;
};

#endif // __anSimulatorEEGSpheres_c_H__
