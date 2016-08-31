//$3    31.07.2001 Matthias D. removed bug in constructor
//$2    24.07.2001 Matthias D. changed members and access functions
//$1    09.08.2000 Matthias D. created

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
// IP_SIM_ParameterSimulator_c.cpp: implementation of the IP_SIM_ParameterSimulator_c class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAnalysisPch.h"

using namespace std;

#include "IP_SIM_ParameterSimulator_c.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

IP_SIM_ParameterSimulator_c::IP_SIM_ParameterSimulator_c(std::string name,std::string unit,double min,double max,double standard, double deviant, double value)
:
m_parameter_name(name),
m_parameter_unit(unit),
m_parameter_minimum(min),
m_parameter_maximum(max),
m_parameter_standard(standard),
m_parameter_deviant(deviant),
m_parameter_value(value),
m_active(true)
{

}

IP_SIM_ParameterSimulator_c::~IP_SIM_ParameterSimulator_c()
{

}

bool IP_SIM_ParameterSimulator_c::setParamterName(std::string &inName)
{
    m_parameter_name = inName;
    return true;
}

std::string IP_SIM_ParameterSimulator_c::getParameterName()
{
    return m_parameter_name;
}

bool IP_SIM_ParameterSimulator_c::setParameterUnits(std::string &inUnits)
{
    m_parameter_unit = inUnits;
    return true;
}

std::string IP_SIM_ParameterSimulator_c::getParameterUnits()
{
    return m_parameter_unit;
}

double IP_SIM_ParameterSimulator_c::getMinimum()
{
    return m_parameter_minimum;
}

bool IP_SIM_ParameterSimulator_c::setMinimum(double inMinimum)
{
    m_parameter_minimum = inMinimum;
    return true;
}

double IP_SIM_ParameterSimulator_c::getMaximum()
{
    return m_parameter_maximum;
}

bool IP_SIM_ParameterSimulator_c::setMaximum(double inMaximum)
{
    m_parameter_maximum = inMaximum;
    return true;
}

double IP_SIM_ParameterSimulator_c::getStandard()
{
    return m_parameter_standard;
}

bool IP_SIM_ParameterSimulator_c::setStandard(double inStandard)
{
    m_parameter_standard = inStandard;
    return true;
}

double IP_SIM_ParameterSimulator_c::getDeviant()
{
    return m_parameter_deviant;
}

bool IP_SIM_ParameterSimulator_c::setDeviant(double inDeviant)
{
    m_parameter_deviant = inDeviant;
    return true;
}

double IP_SIM_ParameterSimulator_c::getValue()
{
    return m_parameter_value;
}

bool IP_SIM_ParameterSimulator_c::setValue(double inValue)
{
    m_parameter_value = inValue;
    return true;
}

bool IP_SIM_ParameterSimulator_c::setValueToStandard()
{
    m_parameter_value = m_parameter_standard;
    return true;
}

bool IP_SIM_ParameterSimulator_c::setValueToDeviant()
{
    m_parameter_value = m_parameter_deviant;
    return true;
}
