//$2    24.02.2000    Frank N.    Changed define for bad_allcoc
//$1    17.08.1999    Frank Z.    Created
//-------------------------------------------------------------------------
//                                      ANT Software BV (c)
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
//-------------------------------------------------------------------------
// File name : CON_UtilitiesDef.h
//-------------------------------------------------------------------------

// CON_UtilitiesDef.h : defintion header file for UtilitiesDll
//

#ifndef __ANALYSISDEF_H__
#define __ANALYSISDEF_H__

#ifndef ANALYSIS_EXPORT

#ifndef _WINDLL

#define ANALYSIS_EXPORT

#else

#ifndef ANALYSIS_DLL
#define ANALYSIS_EXPORT __declspec( dllimport )
#else
#define ANALYSIS_EXPORT __declspec( dllexport )
#endif // ANALYSIS_DLL

#endif // _WINDLL

#endif // ANALYSIS_EXPORT

#endif // __ANALYSISDEF_H__
