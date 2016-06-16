//$8    21.05.2003  Anwander A.  changed setSystemMatrix()
//$7    22.08.2002  Matthias D.  Added setSystemMatrix()
//$6    06.02.2002  Matthias D.  Corrected spelling Error
//$5    03.08.2001  Matthias D.  Added virtual bool DoSetParamterValues() =0;
//$4    20.97.2001  TRK          Some changes with access functions
//$3    18.08.2001  Matthias D.  Added Methods to get and set simulator parameters for error estimation package
//$2    21.08.2000  Matthias D   changed description of input parameters for computeGainMatrix
//$1    01.08.2000  Matthias D   created

#ifndef __IP_SIM_ABSTRACTSIMULATOREEGMEG_C_H__
#define __IP_SIM_ABSTRACTSIMULATOREEGMEG_C_H__

#include "CON_Vector_t.h"
#include "CON_Matrix_t.h"
#include "CON_Block_t.h"


//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche.
//
//    Permission is hereby granted, free of charge, to any person obtaining
//    a copy of the NeuroFEM software and associated documentation files
//    (the "Software"), to deal in the Software without restrictions, including
//    without limitations the rights to use, copy, modify, merge, publish,
//    distribute, sublicense, and/or sell copies of the Software, and to permit
//    persons to whom the Software is furnished to do so, subject to the following
//    conditions:
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//    THE SOFTWARE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
//    SOFTWARE IS WITH YOU.
//
//    The above copyright notice, permission and warranty notices shall be
//    included in all copies of the Software. A copyright notice of the
//    contributing authors and the above permission and warranty notices
//    shall be included in all source code files of the Software.
//
///////////////////////////////////////////////////////////////////////////////
////////////////////// --- abstract base classes for EEG/MEG analysis --- ///////////////////////////////////

// IP_SIM_AbstractSimulatorEEGMEG_c.h: interface for the IP_SIM_AbstractSimulatorEEGMEG_c class.
// Abstract class definition for all kind of EEG/MEG Simulators
//////////////////////////////////////////////////////////////////////


#include "IP_SIM_AbstractSimulator_c.h"
#include "IP_SIM_ParameterSimulator_c.h"

class ANALYSIS_EXPORT IP_SIM_AbstractSimulatorEEGMEG_c : public IP_SIM_AbstractSimulator_c
{
protected:
    IP_SIM_AbstractSimulatorEEGMEG_c();
public:
    virtual ~IP_SIM_AbstractSimulatorEEGMEG_c();

    // A set of methods to set simulator paramters and get parameter values
    int  getNumberofParameters(); // returns number of ACTIVE parameters

    // retrieve parameter by number
    bool ParameterNgetName(    int ParameterNumber, std::string& outName);
    bool ParameterNgetUnit(    int ParameterNumber, std::string& outUnit);
    bool ParameterNgetMinimum( int ParameterNumber, double& outMinimum);
    bool ParameterNgetMaximum( int ParameterNumber, double& outMaximum);
    bool ParameterNgetStandard( int ParameterNumber, double& outStandard);
    bool ParameterNgetDeviant( int ParameterNumber, double& outDeviant);
    bool ParameterNgetValue(   int ParameterNumber, double& outValue);

    // set parameter by number
    bool ParameterNsetStandard( int ParameterNumber, double  inStandard);
    bool ParameterNsetDeviant( int ParameterNumber, double  inDeviant);
    bool ParameterNsetValue(   int ParameterNumber, double  inValue);
    bool ParameterNsetToStandard(int ParameterNumber);
    bool ParameterNsetToDeviant(int ParameterNumber);

    // set parameter by name
    bool ParameterNactivate(std::string& Name);
    bool ParameterNactivateAll();
    bool ParameterNdeactivate(std::string& Name);
    bool ParameterNdeactivateAll();
    bool ParameterNsetStandard( std::string& Name, double  inStandard);
    bool ParameterNsetDeviant( std::string& Name, double  inDeviant);
    bool ParameterNsetValue(   std::string& Name, double  inValue);
    bool ParameterNsetToStandard(std::string& Name);
    bool ParameterNsetToDeviant(std::string& Name);

protected:
    virtual bool DoSetParameterValues() =0; // Transfers parameter values to variables used for calculations

public:
    virtual void computeGainMatrix(int NumberDip, const CON_Matrix_t<double>& inPos, const CON_Matrix_t<double>& inDir, CON_Matrix_t<double>& outSimData) = 0;
    // This function should return a gain matrix, each row referring to a channel and each column
    // to a dipole. If no directions are given, a gain matrix is computed, containing 3 columns
    // for each dipole, representing the standard directions (x,y,z). If only one direction is
    // passed, this one is used for all dipoles.
    //
    // Parameter Description for computeGainMatrix
    // NumberDip    Number of Dipoles
    // inPosition   Matrix with Positions;    number of columns = number of dipoles;
    //                                        number of rows = 3 (x,y,z)
    // inDir        Direction of Dipole;    number of columns = number of dipoles;
    //                                      if number of columns is 1 - all dipoles have same direction
    //                                        if number of rows = 3 (x,y,z) directions in cartesian coordinates,
    //                                        if number of rows =  2, directions in polar coordinates (alpha, theta)
    virtual void computeGainMatrix(const CON_Vector_t<double>& inParameters, CON_Matrix_t<double>& outSimData) = 0;

    virtual void setSystemMatrix(const CON_Matrix_t<double>& SystemMatrix)= 0; // to be override in classes fro simulators, which use a system matrix
    virtual void setSystemMatrix(const CON_Matrix_t<float>& SystemMatrix) = 0;


protected:
    std::vector<IP_SIM_ParameterSimulator_c> m_SimParameters;
};

#endif //__anAbstractSimulatorEEGMEG_c_H__
