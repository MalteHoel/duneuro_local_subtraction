//$7    10.07.2000    Matthias D.     Removed capital letters from file names
//$6    05.07.2000    Frank N.    Adapted for Simbio
//$5    02.09.1999    Frank N.    Took out operators for vector
//$4    23.08.1999    Frank N.    Fixed for BC++
//$3    13.08.1999    Frank N.    Use std::vector
//$2    02.07.1999    Frank N.    Added Serialize()
//$1    02.07.1999    Frank N.    Created
//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche,
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
// CON_AbstractPersistentClass_t.h: interface for the CON_AbstractPersistentClass_t template class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_UTABSTARCTPERSISTENTCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
#define AFX_UTABSTARCTPERSISTENTCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <typeinfo>

#include "CON_AbstractPersistenceArchive_i.h"
#include "CON_PersistentClassRegistry_c.h"

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistentClass_t

template <class T >
class NO_VTABLE CON_AbstractPersistentClass_t
{
public:
    // Serialization
    virtual void Serialize(class CON_AbstractPersistenceArchive_i& ar) = 0;
#ifdef WIN32
    // Operators
    friend class CON_AbstractPersistenceArchive_i& operator<<(class CON_AbstractPersistenceArchive_i& ar, const T* pT);
    friend class CON_AbstractPersistenceArchive_i& operator>>(class CON_AbstractPersistenceArchive_i& ar, T*& pT);
#endif
    // Registration
#ifdef WIN32
        inline static bool RegisterAbstractPersistentClass()
    {
        return CON_PersistentClassRegistry_c::RegisterAbstractPersistentClass( typeid( T ) );
    }
#else
        inline static bool RegisterAbstractPersistentClass()
    {
        ASSERT(FALSE);
        return false;
    }

#endif
};

#endif // !defined(AFX_UTABSTARCTPERSISTENTCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
