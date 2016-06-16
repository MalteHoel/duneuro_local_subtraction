//$9    08.02.2002  Matthias D. Adapted for Linux
//$8    10.07.2000  Matthias D. Removed capital letetrs from file names
//$7    05.07.2000    Frank N.    Adapted for Simbio
//$6    02.09.1999    Frank N.    Took out operators for vector
//$5    23.08.1999    Frank N.    Fixed for BC++
//$4    04.08.1999    Frank N.    Added DECLARE_PERSISTENT_CLASS
//$3    31.07.1999    Frank N.    Added std::
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
// CON_PersistentClass_t.h: interface for the CON_PersistentClass_t template class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CON_PERSISTENTCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
#define AFX_CON_PERSISTENTCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_

#if _MSC_VER > 1000
#pragma once

#endif // _MSC_VER > 1000

#include "CON_AbstractPersistenceArchive_i.h"
#include "CON_PersistentClassRegistry_c.h"

#include <typeinfo>

//////////////////////////////////////////////////////////////////////
// Helper makros

#ifdef WIN32
#define DECLARE_PERSISTENT_CLASS(class_name) \
            friend CON_PersistentClass_t< class_name >;
#else
#define DECLARE_PERSISTENT_CLASS(class_name)
#endif

//////////////////////////////////////////////////////////////////////
// CON_PersistentClass_t

template <class T >
class NO_VTABLE CON_PersistentClass_t
{
public:
    // Serialization
#ifdef WIN32
    virtual void Serialize(class CON_AbstractPersistenceArchive_i& ar) = 0;
#endif
#ifdef WIN32
    // Operators
    friend class CON_AbstractPersistenceArchive_i& operator<<(class CON_AbstractPersistenceArchive_i& ar, const T* pT);
    friend class CON_AbstractPersistenceArchive_i& operator>>(class CON_AbstractPersistenceArchive_i& ar, T*& pT);
#endif
    // Creation
    inline static T* CreateInstance()
    {
        return new T;
    }

    // Registration
#ifdef WIN32
    inline static bool RegisterPersistentClass()
    {
        return CON_PersistentClassRegistry_c::RegisterPersistentClass(typeid(T), (void *(__cdecl *)(void))CON_PersistentClass_t< T >::CreateInstance);
    }
#else
    inline static bool RegisterPersistentClass()
    {
        ASSERT(FALSE);
        return false;
    }
#endif
};

#endif // !defined(AFX_UTPERSISTENTCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
