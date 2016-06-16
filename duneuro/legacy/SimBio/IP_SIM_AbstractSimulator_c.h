//$5    01.08.2000  Matthias D. Removed 1. computeWeightedData, 2. computePenalty, 3. updateParameters, 4. compute progress from Abstract Class Defintion (1.,2. will be defined in abstract class interfaces, 3.4. now in anAnlyserMarqaurdt_c)
//$4    10.07.2000    Matthias D.    Removed Capital letters from file names
//$3    05.07.2000    Frank N.    Adapted for Simbio
//$2    22.11.1999    Frank Z.    made destructor virtual
//$1    26.08.1999    Frank Z.    extended


#ifndef __IP_SIM_ABSTRACTSIMULATOR_C_H__
#define __IP_SIM_ABSTRACTSIMULATOR_C_H__

#include "CON_Vector_t.h"
#include "CON_Matrix_t.h"
#include "CON_Block_t.h"

#include "AnalysisDef.h"

//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche, Dr.Frank Zanow, Frank Neumann.
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

// class IP_SIM_AbstractSimulator_c serves as a base class for all methods used in an inverse analysis method.
// IP_SIM_AbstractSimulator_c makes defines the interfaces needed to implement an inverse analysis method.

class ANALYSIS_EXPORT IP_SIM_AbstractSimulator_c
{
protected:
    IP_SIM_AbstractSimulator_c();
public:
    virtual ~IP_SIM_AbstractSimulator_c();
    virtual bool isValidSimulator() = 0;
    virtual void computeGainMatrix(const CON_Vector_t<double>& inParameters,CON_Matrix_t<double>& outSimData) = 0;
};

#endif //__anAbstractSimulator_c_H__
