//$5    08.02.2002  Matthias D. adapted for Linux
//$5    08.05.2001    Frank N.    Don't use namespace ATL for ATL 7.0
//$4    10.07.2000  Matthias D. Removed captila letters from file names
//$3    05.07.2000    Frank N.    Adapted for Simbio
//$2    03.12.1999    Frank N.    Fixed RegisterPersistentCOMClass()
//$1    02.12.1999    Frank N.    Created
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
// CON_PersistentComClass_t.h: interface for the CON_PersistentComClass_t template class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CON_PERSISTENTCOMCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
#define AFX_CON_PERSISTENTCOMCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "CON_AbstractPersistenceArchive_i.h"
#include "CON_PersistentClassRegistry_c.h"

//////////////////////////////////////////////////////////////////////
// Helper makros

#ifdef WIN32
#define DECLARE_PERSISTENT_COM_CLASS(class_name) \
            friend CON_PersistentComClass_t< class_name >;
#else
#define DECLARE_PERSISTENT_COM_CLASS(class_name)
#endif

//////////////////////////////////////////////////////////////////////
// CON_PersistentComClass_t

template <class T >
class NO_VTABLE CON_PersistentComClass_t
{
public:
    // Serialization
    virtual void Serialize(class CON_AbstractPersistenceArchive_i& ar) = 0;
 #ifdef WIN32
    // Operators
    friend class CON_AbstractPersistenceArchive_i& operator<<(class CON_AbstractPersistenceArchive_i& ar, const T* pT);
    friend class CON_AbstractPersistenceArchive_i& operator>>(class CON_AbstractPersistenceArchive_i& ar, T*& pT);
#endif
    // Creation
#ifdef WIN32
    inline static T* CreateInstance()
    {
        ATL::CComObject<T>* pTCOMObj = NULL;
        VERIFY(SUCCEEDED(ATL::CComObject<T>::CreateInstance(&pTCOMObj)));

        return dynamic_cast<T*>(pTCOMObj);
    }
#endif

    // Registration
#ifdef WIN32
    inline static bool RegisterPersistentCOMClass()
    {
    // First register the normal class name:
        bool bResult = CON_PersistentClassRegistry_c::RegisterPersistentClass(typeid(T), (void *(__cdecl *)(void))CON_PersistentComClass_t< T >::CreateInstance);

    // Second register the class name used by the ATL::CComObject<T>:
        bResult |= CON_PersistentClassRegistry_c::RegisterPersistentClass(typeid( ATL::CComObject< T > ), (void *(__cdecl *)(void))CON_PersistentComClass_t< T >::CreateInstance);

        return bResult;
    }
#else
    inline static bool RegisterPersistentCOMClass()
    {
        ASSERT(FALSE);
        return false;
    }
#endif
};

#endif // !defined(AFX_UTPERSISTENTCOMCLASS_T_H__42067FC5_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
