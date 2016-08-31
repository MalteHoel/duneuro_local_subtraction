//$12    17.10.2003  Anwander A. Changed for gcc 3.3 std::type_info
//$11    08.02.2002  Matthias D. Adapted for stand alone workspace
//$10    22.09.2001    Frank N.    Added REGISTER_PERSISTENT_COM_TEMPLATE_CLASS2
//$9    10.07.2000    Matthias D. Removed capital letters from file names
//$8    05.07.2000  Matthias D. Made some changes to work on Linux
//$7    12.01.2000    Frank N.    Added REGISTER_PERSISTENT_TEMPLATE_CLASS2 and REGISTER_ABSTRACT_PERSISTENT_TEMPLATE_CLASS2
//$6    02.11.1999    Frank N.    Added utPersistentCOMClassRegistrar_t
//$5    22.10.1999    Frank N.    Made m_RegisteredPersistentClasses non static
//$4    02.09.1999    Frank N.    Added REGISTER_PERSISTENT_TEMPLATE_CLASS1 and REGISTER_ABSTRACT_PERSISTENT_TEMPLATE_CLASS1
//$3    05.08.1999    TRK            multmap template got third argument (old STL problem)
//$2    31.07.1999    Frank N.    Added std::
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
// CON_PersistentClassRegistry_c.h: interface for the CON_PersistentClassRegistry_c class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CON_PERSISTENTCLASSREGISTRY_C_H__42067FC6_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
#define AFX_CON_PERSISTENTCLASSREGISTRY_C_H__42067FC6_2F96_11D3_96D4_0080C8FD5716__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <typeinfo>

#include "CON_UtilitiesDef.h"
#include "CON_AbstractPersistentClass_t.h"
#include "CON_PersistentClass_t.h"
#include "CON_PersistentComClass_t.h"

class CON_PersistentClassRegistry_c;


//////////////////////////////////////////////////////////////////////
// Some helpers

UTILITIES_EXPORT CON_PersistentClassRegistry_c* utGetPersistentClassRegistry();
UTILITIES_EXPORT void utSetPersistentClassRegistry(CON_PersistentClassRegistry_c* pNewPersistentClassRegistry);

//////////////////////////////////////////////////////////////////////
// utAbstractPersistentClassRegistrar_t

#ifdef WIN32
#define REGISTER_ABSTRACT_PERSISTENT_CLASS(class_name) \
            utAbstractPersistentClassRegistrar_t< class_name > class_name##AbstractPersistentClassRegistrar;

#define REGISTER_ABSTRACT_PERSISTENT_TEMPLATE_CLASS1(class_name, templateType) \
            utAbstractPersistentClassRegistrar_t< class_name< templateType > > class_name##templateType##AbstractPersistentClassRegistrar;

#define REGISTER_ABSTRACT_PERSISTENT_TEMPLATE_CLASS2(class_name, templateType1, templateType2) \
            utAbstractPersistentClassRegistrar_t< class_name< templateType1 < templateType2 > > > class_name##templateType1##templateType2##AbstractPersistentClassRegistrar;
#endif

template < class T >
class utAbstractPersistentClassRegistrar_t
{
public:
    utAbstractPersistentClassRegistrar_t()
    {
          VERIFY( CON_AbstractPersistentClass_t< T >::RegisterAbstractPersistentClass() );
    }
};

//////////////////////////////////////////////////////////////////////
// utPersistentClassRegistrar_t

#ifdef WIN32
#define REGISTER_PERSISTENT_CLASS(class_name) \
            utPersistentClassRegistrar_t< class_name > class_name##PersistentClassRegistrar;

#define REGISTER_PERSISTENT_TEMPLATE_CLASS1(class_name, templateType) \
            utPersistentClassRegistrar_t< class_name< templateType > > class_name##templateType##PersistentClassRegistrar;

#define REGISTER_PERSISTENT_TEMPLATE_CLASS2(class_name, templateType1, templateType2) \
            utPersistentClassRegistrar_t< class_name< templateType1 < templateType2 > > > class_name##templateType1##templateType2##PersistentClassRegistrar;
#endif

template < class T >
class utPersistentClassRegistrar_t
{
public:
    utPersistentClassRegistrar_t()
    {
          VERIFY( CON_PersistentClass_t< T >::RegisterPersistentClass() );
    }
};

//////////////////////////////////////////////////////////////////////
// utPersistentCOMClassRegistrar_t

#ifdef WIN32
#define REGISTER_PERSISTENT_COM_CLASS(class_name) \
            utPersistentCOMClassRegistrar_t< class_name > class_name##PersistentCOMClassRegistrar;

#define REGISTER_PERSISTENT_COM_TEMPLATE_CLASS1(class_name, templateType) \
            utPersistentCOMClassRegistrar_t< class_name< templateType > > class_name##templateType##PersistentCOMClassRegistrar;

#define REGISTER_PERSISTENT_COM_TEMPLATE_CLASS2(class_name, templateType1, templateType2) \
            utPersistentCOMClassRegistrar_t< class_name< templateType1 < templateType2 > > > class_name##templateType1##templateType2##PersistentCOMClassRegistrar;
#endif

template < class T >
class utPersistentCOMClassRegistrar_t
{
public:
    utPersistentCOMClassRegistrar_t()
    {
        VERIFY( CON_PersistentComClass_t< T >::RegisterPersistentCOMClass() );
    }
};

//////////////////////////////////////////////////////////////////////
// CON_PersistentClassRegistry_c

class UTILITIES_EXPORT CON_PersistentClassRegistry_c
{
public:
    CON_PersistentClassRegistry_c();
    virtual ~CON_PersistentClassRegistry_c();

    // Register
#ifdef WIN32
    static bool RegisterPersistentClass(const type_info& typeInfo, void*(__cdecl* pCreateInstanceFunction)());
#else
    static bool RegisterPersistentClass(const std::type_info& typeInfo, void*(pCreateInstanceFunction)());
#endif
    static bool RegisterAbstractPersistentClass(const std::type_info& typeInfo);
#ifdef WIN32
    // Create an instance from the given type
    static void* CreateInstance(const type_info& typeInfo);
    static void* CreateInstance(const std::string& typeName);
#endif
// Data
protected:
    std::multimap<std::string, void*,std::less<std::string> >    m_RegisteredPersistentClasses;
};

#endif // !defined(AFX_UTPERSISTENTCLASSREGISTRY_C_H__42067FC6_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
