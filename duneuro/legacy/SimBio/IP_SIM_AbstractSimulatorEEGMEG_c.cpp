//$4 06.02.2002 Matthias D.  Corrected spelling Error
//$3 03.08.2001 Matthias D.  Added DoSetParamterValues to methods setting simulator parameters
//$2 18.08.2001 Matthias D.  Added Methods to get and set simulator parameters for error estimation package
//$1 01.08.2000 Matthias D.  created

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
// IP_SIM_AbstractSimulatorEEGMEG_c.cpp: implementation of the IP_SIM_AbstractSimulatorEEGMEG_c class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAnalysisPch.h"

using namespace std;

#include "IP_SIM_AbstractSimulatorEEGMEG_c.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

IP_SIM_AbstractSimulatorEEGMEG_c::IP_SIM_AbstractSimulatorEEGMEG_c():
IP_SIM_AbstractSimulator_c()
{

}

IP_SIM_AbstractSimulatorEEGMEG_c::~IP_SIM_AbstractSimulatorEEGMEG_c()
{

}

///////////////////////////////////////////////////////////////////////
// Methods to get and set simulator parameters
///////////////////////////////////////////////////////////////////////

int IP_SIM_AbstractSimulatorEEGMEG_c::getNumberofParameters() // returns number of ACTIVE parameters
{
    int nNumberOfActiveParameters = 0;

    for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].isActive()) nNumberOfActiveParameters++;
        }

    return nNumberOfActiveParameters;
}

// retrieve parameter by number

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetName(int ParameterNumber, std::string& outName)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outName=m_SimParameters[ParameterNumber].getParameterName();
    return true;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetUnit(int ParameterNumber, std::string& outUnit)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outUnit=m_SimParameters[ParameterNumber].getParameterUnits();
    return true;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetMinimum(int ParameterNumber, double& outMinimum)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outMinimum = m_SimParameters[ParameterNumber].getMinimum();
    return true;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetMaximum(int ParameterNumber, double& outMaximum)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outMaximum = m_SimParameters[ParameterNumber].getMaximum();
    return true;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetStandard(int ParameterNumber, double& outStandard)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outStandard = m_SimParameters[ParameterNumber].getStandard();
    return true;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetDeviant(int ParameterNumber, double& outDeviant)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outDeviant = m_SimParameters[ParameterNumber].getDeviant();
    return true;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNgetValue(int ParameterNumber, double& outValue)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    outValue = m_SimParameters[ParameterNumber].getValue();
    return true;
}


// set parameter by number

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetStandard( int ParameterNumber, double  inStandard)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    if (!m_SimParameters[ParameterNumber].setStandard(inStandard)) return false;
    return true;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetDeviant( int ParameterNumber, double  inDeviant)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    if (!m_SimParameters[ParameterNumber].setDeviant(inDeviant)) return false;
    return true;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetValue(   int ParameterNumber, double  inValue)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    if (!m_SimParameters[ParameterNumber].setValue(inValue)) return false;

    if (! DoSetParameterValues()) return false;

    return true;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetToStandard(int ParameterNumber)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    if (!m_SimParameters[ParameterNumber].setValueToStandard()) return false;

    if (! DoSetParameterValues()) return false;

    return true;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetToDeviant(int ParameterNumber)
{
    if ((ParameterNumber < 0)|| (ParameterNumber >= m_SimParameters.size())) return false;

    if (!m_SimParameters[ParameterNumber].setValueToDeviant()) return false;

    if (! DoSetParameterValues()) return false;

    return true;
}


// set parameter by name

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNactivate(std::string& Name)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             m_SimParameters[i].activate();
             return true;
            }
        }
    return false;
}



bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNdeactivate(std::string& Name)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             m_SimParameters[i].deactivate();
             return true;
            }
        }
    return false;
}


bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetStandard( std::string& Name, double  inStandard)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             if(m_SimParameters[i].setStandard(inStandard)) return true;
            }
        }
    return false;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetDeviant( std::string& Name, double  inDeviant)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             if(m_SimParameters[i].setDeviant(inDeviant)) return true;
            }
        }
    return false;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetValue(   std::string& Name, double  inValue)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             if(m_SimParameters[i].setValue(inValue))
                {
                 if (! DoSetParameterValues()) return false;
                 return true;
                }
            }
        }
    return false;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetToStandard(std::string& Name)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             if(m_SimParameters[i].setValueToStandard()) return true;
            }
        }
    return false;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNsetToDeviant(std::string& Name)
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
         if ( m_SimParameters[i].getParameterName() == Name )
            {
             if(m_SimParameters[i].setValueToDeviant()) return true;
            }
        }
    return false;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNactivateAll()
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
           m_SimParameters[i].activate();
        }

    return true;
}

bool IP_SIM_AbstractSimulatorEEGMEG_c::ParameterNdeactivateAll()
{
     for(int i =0; i< m_SimParameters.size();i++)
        {
           m_SimParameters[i].deactivate();
        }

    return true;
}
