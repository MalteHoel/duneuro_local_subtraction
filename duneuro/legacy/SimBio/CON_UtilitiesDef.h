//$19   02.01.2006  Graichen U. removed bugs for new gnu compiler
//$18    13.02.2003  Anwander A. Changer sqrtl to sqrt for portability
//$17   29.04.2002  Stefan      Defined a macro for memory leaks localisation
//$16    08.02.2002  Matthias D. Changed for Linux
//$15   23.06.2000  Matthias D  Changes for ASSERT and VERIFY derfinitions
//$14   15.06.2000  Matthias D  Added some typdefs for SDK and MFC types to make it work on Linux
//$13    24.02.2000    Frank N.    Changed define for bad_allcoc
//$12    03.09.1999    Frank N.    Suppress more STL related warnings for BC++
//$11    23.08.1999    Frank N.    Fixed for BC++
//$10    22.08.1999    Frank N.    Suppress a couple of STL specific warnings
//$9    05.08.1999    TRK            BC++ redefine for xalloc
//$8    06.08.1999    Frank Z.    Made changes working under MSC
//$8    05.08.1999    TRK            Einige defines verndert, damit es unter Borland kompiliert
//$7    31.07.1999    Frank N.    Added #include <set>
//$6    14.07.1999    Frank N.    Added #include <shlwapi.h>
//$5    05.07.1999    Frank N.    Added #include <objbase.h>
//$4    01.07.1999    Frank N.    Added NO_VTABLE and #include <map>                // STL map class
//$3    14.05.1999    Frank N.    added <vector> include
//$2    10.05.1999    Frank N.    changed includes for STL
//$1    05.05.1999    Frank N.    completly redone
//-------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche, Dr. U. Graichen
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

#ifndef _CON_UTILITIESDEF_H_
#define _CON_UTILITIESDEF_H_

#ifdef _MSC_VER
// disable warning C4786: symbol greater than 255 character,// okay to ignore
#pragma warning(disable: 4786)

#ifndef MAXPATH
#define MAXPATH _MAX_PATH
#endif

#endif    // _MSC_VER


#define MAX_LABEL_LENGTH 100
//aa 13.02.2002 sqrtl is optional in ANSI C and don't works on SGI !!!
//#define sqrx(x) ((x<=1e-36)?0:sqrtl(x)) // safe square root
#define sqrx(x) ((x<=1e-36)?0:sqrt(x)) // safe square root

// Some typedefs that make it work on Linux Matthias D. 15.06.2000
#ifdef WIN32

#define NO_VTABLE __declspec(novtable)

#else

#define NO_VTABLE

#include <assert.h>
#include <string.h>
#include <string>

typedef unsigned int UINT;
typedef unsigned char BYTE;
typedef unsigned short WORD;
typedef long int LONG;
typedef unsigned int DWORD;

#define VERIFY(f)   assert(f)
#define ASSERT(f)   assert(f)

#define ZeroMemory(Destination,Length) memset((Destination),0,(Length))

#ifndef FALSE
#define FALSE               0
#endif

#ifndef TRUE
#define TRUE                1
#endif

#define _T(x)               x

/* unreachable code, #ifdef WIN32 && #ifndef WIN32 */
/* #ifdef WIN32 */
/* class runtime_error: public exception{ */
/* public: */
/*     explicit runtime_error(const string& _S) */
/*         : exception(), _Str(_S) {}; */
/*     virtual ~runtime_error() */
/*                 {} */
/*         virtual const char *what() const */
/*                 {return (_Str.c_str()); } */

/* protected: */
/*        // virtual void _Doarise() const */
/*       //        {_RAISE(*this):} */
/* private: */
/*          string _Str; */
/*  }; */
/* #endif */

#endif

#ifndef UTILITIES_EXPORT

#ifndef _WINDLL

#define UTILITIES_EXPORT

#else

#ifndef UTILITIES_DLL
#define UTILITIES_EXPORT __declspec( dllimport )
#else
#define UTILITIES_EXPORT __declspec( dllexport )
#endif // UTILITIES_DLL

#endif // _WINDLL

#endif // UTILITIES_EXPORT

#endif // _UTILITIESDEF_H_

//*****************************************************************************************
// helper for memory leaks localisation
//
// the macro 'TRACE_ALLOC" put a trace on the output view of Visual Studio for
// mem-allocation of an object.
// This macro is default disabled, because it makes the program very slow in debug mode
//
// If you want to use the macro, enable the define 'MEM_ANALYZER' on your local machine,
// without checking out this file. After this, a rebuild is required!
//*****************************************************************************************

//#define MEM_ANALYZER

#ifdef WIN32
    #ifdef MEM_ANALYZER
        #define TRACE_ALLOC(pointer, func_name) TRACE3(_T("%d bytes at 0x%x in %s allocated\n"), sizeof(*(pointer)), (pointer), (func_name))
    #else
        #define TRACE_ALLOC(pointer, func_name) // do nothing :-)
    #endif // MEM_ANALYZER
#else
    #define TRACE_ALLOC(pointer, func_name) // do nothing :-)
#endif // WIN32
