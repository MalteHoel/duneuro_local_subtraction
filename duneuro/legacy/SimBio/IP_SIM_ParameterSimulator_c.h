//$2    20.07.2001 TRK           introduced standard and deviant, bugs in argument lists, constructor
//$1    09.08.2000 Matthias D. created


#ifndef __IP_SIM_PARAMETERSIMULATOR_C_H__
#define __IP_SIM_PARAMETERSIMULATOR_C_H__

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
// IP_SIM_ParameterSimulator_c.h: interface for the IP_SIM_ParameterSimulator_c class.
//
//////////////////////////////////////////////////////////////////////



class ANALYSIS_EXPORT IP_SIM_ParameterSimulator_c
{
public:
    IP_SIM_ParameterSimulator_c(std::string name,std::string unit,double min,double max,double standard, double deviant, double value);
    virtual ~IP_SIM_ParameterSimulator_c();

    bool setParamterName(std::string &inName);
    std::string getParameterName();
    bool setParameterUnits(std::string &inUnits);
    std::string getParameterUnits();
    double getMinimum();
    bool   setMinimum(double inMinimum) ;
    double getMaximum();
    bool   setMaximum(double inMaximum) ;
    double getStandard();
    bool   setStandard(double inStandard);
    double getDeviant();
    bool   setDeviant(double inDeviant);
    double getValue();
    bool setValue(double inValue);
    bool setValueToStandard();
    bool setValueToDeviant();
    void activate() {m_active = true;}
    void deactivate() {m_active = false;}
    bool isActive() {return m_active;}

protected:
std::string        m_parameter_name;
std::string        m_parameter_unit;
double            m_parameter_minimum;
double            m_parameter_maximum;
double            m_parameter_standard;
double            m_parameter_deviant;
double            m_parameter_value;
bool            m_active;
};

#endif // __anParameterSimulator_c_H__
