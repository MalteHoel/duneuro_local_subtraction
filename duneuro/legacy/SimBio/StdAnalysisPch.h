//$14   29.04.2002  Stefan      moved '#define TRACE_ALLOC' to a global place --> CON_UTILITIESDEF.h
//$13    26.04.2002    Velko        added #define TRACE_ALLOC
//$12    07.02.2002  Matthias D. Adapted for Linux
//$11    18.01.2002  Matthias D. Added #include <SeizureDetectorDef.h>
//$10    08.05.2001    Frank N.    Don't use namespace ATL for ATL 7.0
//$9    10.07.2000  Matthias D. Removed Capital letters from file names
//$8    05.07.2000    Frank N.    Adapted for Simbio
//$7    16.06.2000    Frank N.    Added #include <SpikeDetectorDef.h>
//$6    09.06.2000    Frank N.    Took out include for data core
//$5    28.03.2000    Jose F.        Added more STL includes ( functional  & algorithm )
//$4    28.03.2000    Ralf        Added afxcmn.h
//$3    18.03.2000    Frank N.    Added include "atlbase.h"
//$2    15.03.2000    Frank N.    Added more STL includes
//$1    22.12.1999    Frank N.    Created
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
// StdAnalysisPch.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDANALYSIS_H__C6247D9D_BFA2_11D2_9B60_3295D3000000__INCLUDED_)
#define AFX_STDANALYSIS_H__C6247D9D_BFA2_11D2_9B60_3295D3000000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef WIN32

#include <valarray>            // STL valarray class

#define VC_EXTRALEAN        // Exclude rarely-used stuff from Windows headers
#define _ATL_NO_AUTOMATIC_NAMESPACE // Don't use automatically the ATL namespace (ATL 7.0)

#include <afxwin.h>         // MFC core and standard components
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>            // MFC support for Windows Common Controls
#endif // _AFX_NO_AFXCMN_SUPPORT


#include <objbase.h>
#include <shlwapi.h>

#endif

// disable warning C4786: symbol greater than 255 character,// okay to ignore
#pragma warning(disable: 4786)

#include <string>            // STL string class
#include <vector>            // STL vector
#include <map>                // STL map class
#include <set>                // STL set class
#include <list>                // STL list class
#include <functional>
#include <algorithm>
#include <exception>

// Utilities library components
#include "CON_UtilitiesDef.h"

#ifdef WIN32
// SpikeDetector library components
#include <spikedetectordef.h>

// SpikeDetector library components
#include <seizuredetectordef.h>
#endif // WIN32

// Analysis library components
#include "AnalysisDef.h"


#ifdef WIN32

// ATL
#include "atlbase.h"

#endif // WIN32


//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDANALYSIS_H__C6247D9D_BFA2_11D2_9B60_3295D3000000__INCLUDED_)
