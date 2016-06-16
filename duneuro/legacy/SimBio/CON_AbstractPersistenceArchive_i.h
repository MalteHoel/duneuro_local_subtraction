//$38   02.01.2006  Graichen U. removed bugs for new gnu compiler (added ifdef WIN32)
//$37    17.10.2003  Anwander A. adapted for gcc 3.3 added typename
//$36   19.03.2002  Matthias D. Adapted for Linux
//$35   08.02.2002  Matthias D. Adapted for Linux
//$34    02.11.2001  Maurice        Added member m_bToClipboard, added access functions and changed constructors
//$33    17.10.2001    Ben            Added ASSERT() to all utLoadXXXX()-style functions
//$32    03.08.2001    Ben            Fixed memory leaks in utLoadVectorOfObjectReferences() and in utLoadMapWithObjectsOfVectorWithObjects()
//$31    11.07.2001    Ben            Added new collector operator functions
//$30    08.05.2001    Frank N.    Don't use 'using namespace ATL;' for ATL 7.0
//$29    29.01.2001    Frank N.    Added GetLength() and GetPosition()
//$28    11.12.2000    Frank N.    Added utStore/LoadSet()
//$27    12.10.2000  Matthias D. Fixed another problem with an ifndef
//$26    03.08.2000    Frank N.    Fixed another problem with an ifndef
//$25   10.07.2000  Matthias D. Removed capital letters from file names
//$24   05.07.2000  Matthias D. Added inculde-file typeinfo
//$23    24.03.2000    Ben            Inserted to ASSERTS( m.find(key) == m.end() )
//$22    18.03.2000    Frank N.    Implemented operator<</>> for CComVariant
//$21    24.02.2000    Frank N.    Changed define for bad_allcoc
//$20    15.02.2000    Frank N.    Changed return type of GetVersionNumber()
//$19    12.01.2000    Frank Z.    Added utStoreList/utLoadList
//$18    12.01.2000    Frank N.    Added utStore/LoadMapOfObjectPointers()
//$17    23.12.1999    Frank N.    Replaced _MSC_VER by _ASA
//$16    26.10.1999    Frank N.    Implemented operator<</>> for GUID
//$15    03.09.1999    Frank N.    Don't use the operator<</>> for stl collectors anymore
//$14    02.09.1999    Frank N.    Map object pointers
//$13    20.08.1999    Frank N.    Fixed operator<< for map
//$12    18.08.1999    Frank N.    1. Added streaming for WINDOWS primitives 2. Added specialized operators for vector<string>
//$11    19.08.1999    TRK            added stl versions without default arguments
//$10   17.08.1999    Frank N.    Added operators for list
//$9    17.08.1999    Frank N.    Added operators for map
//$8    10.08.1999    Frank N.    << operator for bool global existent
//$7    06.08.1999    Frank Z.    << operator for bool compiler-dependent existent
//$6    05.08.1999    TRK            << operator fr bool wieder eingefhrt
//$5    05.08.1999    TRK            bei einigen Funktionen inline statement entfernt (Borland wollte es so)
//$4    31.07.1999    Frank N.    Removed using namespace std directive (it's very bad to have it in a header file)
//$3    29.07.1999    Ben            Added "using namespace std" directive
//$2    02.07.1999    Frank N.    Added header information
//$1    02.07.1999    Frank N.    Created
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
/////////////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i.h: interface for the CON_AbstractPersistenceArchive_i class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CON_ABSTRACTPERSISTENCEARCHIVE_I_H__42067FC7_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
#define AFX_CON_ABSTRACTPERSISTENCEARCHIVE_I_H__42067FC7_2F96_11D3_96D4_0080C8FD5716__INCLUDED_

#include <stdexcept>

#include <string>
#include <vector>
#include <map>
#include <list>
#include <set>

#include <typeinfo>

#include <time.h>

#include "CON_UtilitiesDef.h"
#include "CON_PersistentClassRegistry_c.h"


// BL: Replacement of BOOL by bool had to be done throughout SimBio
//     in order to get MPICH2 running with SimBio.
//#ifndef WIN32
////#define  BOOL bool
//#endif

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CON_AbstractPersistenceArchive_i;

template <class T>
class utSmartPtr_t;

/////////////////////////////////////////////////////////////////////////////
// Compression modes of the archive

enum utCompressionMode_e
{
    CM_UNCOMPRESSED        = 0    // no compression used
};

/////////////////////////////////////////////////////////////////////////////
// Encrytion modes of the archive

enum utEncryptionMode_e
{
    EM_NOTENCRYPTED        = 0    // no encryption used
};

////////////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i support for polymorphic reading/writing of objects

// Pointer mapping constants
#define dwNullTag        ((DWORD)0)                // special tag indicating NULL ptrs
#define dwClassNameTag    ((DWORD)0xFFFFFFFF)        // special tag indicating class definition
#define nMaxMapCount    ((DWORD)0x3FFFFFFE)        // 0x3FFFFFFE last valid mapCount

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i Template operator functions

class UTILITIES_EXPORT CON_AbstractPersistenceArchive_i
{
protected:
    CON_AbstractPersistenceArchive_i(bool bLoading, int nVersion = -1, utCompressionMode_e compressionMode = CM_UNCOMPRESSED, utEncryptionMode_e encryptionMode = EM_NOTENCRYPTED, bool bToClipboard = FALSE);

public:
    virtual ~CON_AbstractPersistenceArchive_i();

public:
// Attributes
    inline bool IsLoading() const { return m_bLoading; }
    inline bool IsStoring() const { return !IsLoading(); }

    inline bool toClipboard() const {return m_bToClipboard;}
    inline void setToClipboard(bool bToClipboard){m_bToClipboard = bToClipboard;}

    time_t GetCreationTime() const { return m_ctime; }
    time_t GetLastWriteTime() const { return m_mtime; }
    time_t GetLastAccessTime() const { return m_atime; }
    inline int GetVersionNumber() const { return m_nVersionNumber; }
    inline utCompressionMode_e GetCompressionMode() const { return m_eCompressionMode; }
    inline utEncryptionMode_e GetEncryptionMode() const { return m_eEncryptionMode; }

    virtual bool ContainsOnlyHeader() const = 0;
    DWORD GetHeaderLength() const;
    virtual DWORD GetLength() const = 0;
    virtual DWORD GetPosition() const = 0;

public:
// Operations:
    virtual UINT Read(void* lpBuf, UINT nMax) = 0;
    virtual void Write(const void* lpBuf, UINT nMax) = 0;
    virtual void Close() = 0;
    virtual void Abort() = 0;   // close and shutdown without exceptions

    // insertion operations
    virtual CON_AbstractPersistenceArchive_i& operator<<(bool b) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(BYTE by) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(WORD w) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(LONG l) = 0;
#ifdef WIN32
    virtual CON_AbstractPersistenceArchive_i& operator<<(DWORD dw) = 0;
#endif
    virtual CON_AbstractPersistenceArchive_i& operator<<(float f) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(double d) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(int i) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(short w) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(char ch) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(unsigned u) = 0;

    virtual CON_AbstractPersistenceArchive_i& operator<<(const std::string& str) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator<<(const struct tm& tm) = 0;
#ifdef WIN32
    virtual CON_AbstractPersistenceArchive_i& operator<<(const ATL::CComVariant& var) = 0;
#endif

    // extraction operations
    virtual CON_AbstractPersistenceArchive_i& operator>>(bool& b) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(BYTE& by) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(WORD& w) = 0;
#ifdef WIN32
    virtual CON_AbstractPersistenceArchive_i& operator>>(DWORD& dw) = 0;
#endif
    virtual CON_AbstractPersistenceArchive_i& operator>>(LONG& l) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(float& f) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(double& d) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(int& i) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(short& w) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(char& ch) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(unsigned& u) = 0;

    virtual CON_AbstractPersistenceArchive_i& operator>>(std::string& str) = 0;
    virtual CON_AbstractPersistenceArchive_i& operator>>(struct tm& tm) = 0;
#ifdef WIN32
    virtual CON_AbstractPersistenceArchive_i& operator>>(ATL::CComVariant& var) = 0;
#endif

public:
// Access to the object maps
    long FindInStoreMap(void* pObj) const;
    bool AddToStoreMap(void* pObj);
#ifdef WIN32
    void* FindInLoadArray(long nIndex) const;
#endif
    bool AddToLoadArray(void* pObj);

protected:
// Object mapping to ensure uniqueness:
    void CheckObjectCount();  // throw exception if m_nMapCount is too large
    bool InitializeObjectMaps();

protected:
// Header:
    bool WriteHeader();
    bool ReadHeader();

protected:
// Header information (stored):
    // General information:
    time_t                            m_ctime;          // creation date/time of file
    time_t                            m_mtime;          // last modification date/time of file
    time_t                            m_atime;          // last access date/time of file
    // Versioning for the whole archive:
    int                                m_nVersionNumber;
    // Compression
    utCompressionMode_e                m_eCompressionMode;
    // Encryption
    utEncryptionMode_e                m_eEncryptionMode;

// Run-time only:
    // Loading/Storing
    bool                            m_bLoading;
    bool                            m_bToClipboard;

    // Maps for mapping between:
    UINT                            m_nMapCount;
    union
    {
    // 1. Pointer in stream to application object
        std::vector<void*>*            m_pLoadArray;
    // 2. Pointer of application object to pointer in stream
        std::map<void*, DWORD>*        m_pStoreMap;
    };
};

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i Template operator functions

template< class T >
CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const T* pObj)
{
// Write this object to the archive:
    long nObjIndex = 0;

// Make sure we are storing:
    ASSERT( ar.IsStoring() );

// Make sure that we are on the right machine:
    ASSERT(sizeof(dwNullTag) == 4);
    ASSERT(sizeof(dwClassNameTag) == 4);

    if (pObj == NULL)
    {
    // Save out null tag to represent NULL pointer
        ar << dwNullTag;
    }
    else if ( (nObjIndex = ar.FindInStoreMap( const_cast<T*>(pObj) )) != -1 )
    {
    // Save out index of already stored object
        ar << nObjIndex;
    }
    else
    {
    // Write class of object first
        const std::type_info& typeInfo = typeid(*pObj);
        std::string strTypeName = typeInfo.name();
        ar << dwClassNameTag;
        ar << strTypeName;

    // Enter in stored object table, checking for overflow
        bool bResult = ar.AddToStoreMap( const_cast<T*>(pObj) );
        ASSERT(bResult);

        if (!bResult)
        {
        // An error occured
#ifdef _AFXDLL
            AfxThrowNotSupportedException();
#else
            throw std::runtime_error(_T("notSupported"));
#endif
        }

    // Cause the object to serialize itself
        const_cast<T*>(pObj)->Serialize(ar);
    }

    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, T*& pObj)
{
#ifdef WIN32
// Read this object from the archive:
    void* theObj = NULL;
    ASSERT(pObj == NULL);

// Make sure we are loading:
    ASSERT( ar.IsLoading() );

    ASSERT(dwNullTag == 0);

    DWORD theTag;
    ar >> theTag;
    if (theTag == dwNullTag)
    {
        pObj = NULL;
    }
    else if ( (theObj = ar.FindInLoadArray(theTag)) != NULL )
    {
#ifdef _DEBUG
        pObj = reinterpret_cast<T*>(theObj);
        ASSERT(pObj != NULL);
        if (pObj == NULL)
        {
#ifdef _AFXDLL
          AfxThrowArchiveException(CArchiveException::badClass, _T(""));
#else
            throw std::runtime_error(_T("badClass"));
#endif
        }
#else
        pObj = static_cast<T*>(theObj);
#endif
    }
    else if (theTag == dwClassNameTag)
    {
        std::string strTypeName;

        ar >> strTypeName;

    // Create the appropriate object:
        pObj = static_cast< T* >( CON_PersistentClassRegistry_c::CreateInstance(strTypeName) );
        ASSERT(pObj != NULL);

        if (pObj == NULL)
        {
#ifdef _AFXDLL
            AfxThrowMemoryException();
#else
            throw new std::bad_alloc;
#endif
        }

    // Add to mapping array BEFORE de-serializing
        if ( !ar.AddToLoadArray(pObj) )
        {
#ifdef _AFXDLL
            AfxThrowArchiveException(CArchiveException::badIndex, _T(""));
#else
            throw std::runtime_error(_T("badIndex"));
#endif
        }

    // Serialize the object
        pObj->Serialize(ar);
    }
    else
    {
    // This should never happen:
        ASSERT(false);
#ifdef _AFXDLL
        AfxThrowArchiveException(CArchiveException::badSchema, _T(""));
#else
            throw std::runtime_error(_T("badSchema"));
#endif
    }
#endif //WIN32
    return ar;
}

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i smart pointer operator functions

template< class T >
CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const utSmartPtr_t< T >& pSmartPtr)
{
    ar << (T*)pSmartPtr;

    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, utSmartPtr_t< T >& pSmartPtr)
{
    T* pObj = NULL;
    ar >> pObj;
    pSmartPtr = pObj;

    return ar;
}

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i std::pair operator functions

template< class T1, class T2 >
CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const std::pair< T1, T2 >& Pair)
{
    ar << Pair.first;
    ar << Pair.second;

    return ar;
}

template< class T1, class T2 >
CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, std::pair< T1, T2 >& Pair)
{
    ar >> Pair.first;
    ar >> Pair.second;

    return ar;
}

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i Templatized collector operator functions

template< class T >
CON_AbstractPersistenceArchive_i& utStoreVectorOfObjectPointers(CON_AbstractPersistenceArchive_i& ar, const std::vector< T >& v)
{
#ifdef WIN32
    ar << v.size();

    for (typename std::vector< T >::const_iterator theIterator = v.begin(); theIterator != v.end(); theIterator++)
    {
        ar << *theIterator;
    }
#endif
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utLoadVectorOfObjectPointers(CON_AbstractPersistenceArchive_i& ar, std::vector< T >& v)
{
#ifdef WIN32
    // make sure to pass an empty vector
    ASSERT( v.empty() );


    typename std::vector< T >::size_type nSize = 0;
    ar >> nSize;

    v.reserve(nSize);

    for (typename std::vector< T >::size_type i = 0; i < nSize; i++)
    {
        typename std::vector< T >::value_type CurrElement = NULL;
        ar >> CurrElement;

        v.push_back(CurrElement);
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utStoreVectorOfObjectReferences(CON_AbstractPersistenceArchive_i& ar, const std::vector< T >& v)
{
#ifdef WIN32
    ar << v.size();

    for (typename std::vector< T >::size_type i = 0; i < v.size(); i++)
    {
        typename std::vector< T >::const_reference currElement = v.at(i);
        ar << &currElement;
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utLoadVectorOfObjectReferences(CON_AbstractPersistenceArchive_i& ar, std::vector< T >& v)
{
#ifdef WIN32
    // make sure to pass an empty vector
    ASSERT( v.empty() );

    typename std::vector< T >::size_type nSize = 0;
    ar >> nSize;

    v.reserve(nSize);

    for (typename std::vector< T >::size_type i = 0; i < nSize; i++)
    {
        typename std::vector< T >::value_type* pCurrElement = NULL;
        ar >> pCurrElement;
        ASSERT(pCurrElement != NULL);

        if (pCurrElement != NULL)
        {
            v.push_back(*pCurrElement);

            // push_back makes a copy of currentElement
            // delete object
            delete pCurrElement;
        }
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utStoreVector(CON_AbstractPersistenceArchive_i& ar, const std::vector< T >& v)
{
#ifdef WIN32
    ar << v.size();

    for (typename std::vector< T >::size_type i = 0; i < v.size(); i++)
    {
        typename std::vector< T >::const_reference currElement = v.at(i);
        ar << currElement;
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utLoadVector(CON_AbstractPersistenceArchive_i& ar, std::vector< T >& v)
{
#ifdef WIN32
    // make sure to pass an empty vector
    ASSERT( v.empty() );


    typename std::vector< T >::size_type nSize = 0;
    ar >> nSize;

    v.reserve(nSize);

    for (typename std::vector< T >::size_type i = 0; i < nSize; i++)
    {
        typename std::vector< T >::value_type CurrElement;
        ar >> CurrElement;

        v.push_back(CurrElement);
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utStoreList(CON_AbstractPersistenceArchive_i& ar, const std::list< T >& l)
{
#ifdef WIN32
    ar << l.size();

    for (typename std::list< T >::const_iterator theIterator = l.begin(); theIterator != l.end(); theIterator++)
    {
        typename std::list< T >::const_reference currElement = *theIterator;
        ar << currElement;
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utLoadList(CON_AbstractPersistenceArchive_i& ar, std::list< T >& l)
{
#ifdef WIN32
    // make sure to pass an empty list
    ASSERT( l.empty() );

    typename std::list< T >::size_type nSize = 0;
    ar >> nSize;

    for (typename std::list< T >::size_type i = 0; i < nSize; i++)
    {
        typename std::list< T >::value_type CurrElement;
        ar >> CurrElement;

        l.push_back(CurrElement);
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utStoreListOfObjectPointers(CON_AbstractPersistenceArchive_i& ar, const std::list< T >& l)
{
#ifdef WIN32
    ar << l.size();

    for (typename std::list< T >::const_iterator theIterator = l.begin(); theIterator != l.end(); theIterator++)
    {
        ar << *theIterator;
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utLoadListOfObjectPointers(CON_AbstractPersistenceArchive_i& ar, std::list< T >& l)
{
#ifdef WIN32
    // make sure to pass an empty list
    ASSERT( l.empty() );


    typename std::list< T >::size_type nSize = 0;
    ar >> nSize;

    for (typename std::list< T >::size_type i = 0; i < nSize; i++)
    {
        typename std::list< T >::value_type CurrElement;
        ar >> CurrElement;

        l.push_back(CurrElement);
    }
#endif //WIN32
    return ar;
}

template< class Key, class T >
CON_AbstractPersistenceArchive_i& utStoreMap(CON_AbstractPersistenceArchive_i& ar, const std::map< Key, T >& m)
{
#ifdef WIN32
    ar << m.size();

    for (typename std::map< Key, T >::const_iterator theIterator = m.begin(); theIterator != m.end(); ++theIterator)
    {
        ar << theIterator->first;
        ar << theIterator->second;
    }
#endif //WIN32
    return ar;
}

template< class Key, class T >
CON_AbstractPersistenceArchive_i& utLoadMap(CON_AbstractPersistenceArchive_i& ar, std::map< Key, T >& m)
{
#ifdef WIN32
    // make sure to pass an empty map
    ASSERT( m.empty() );

    typename std::map< Key, T >::size_type nSize = 0;
    ar >> nSize;

    for (typename std::map< Key, T >::size_type i = 0; i < nSize; i++)
    {
        typename std::map< Key, T >::key_type key;
#ifdef WIN32
        std::map< Key, T >::referent_type referent;
#else
        typename std::map< Key, T >::reference referent;
#endif
        ar >> key;
        ar >> referent;

        ASSERT( m.find( key ) == m.end() );
#ifdef WIN32
        m.insert( std::pair<std::map< Key, T >::key_type, std::map< Key, T >::referent_type>(key, referent));
#else
           m.insert( typename std::pair<typename std::map< Key, T >::key_type, typename std::map< Key, T >::reference>(key, referent));
#endif
    }
#endif //WIN32
    return ar;
}

template< class Key, class T >
CON_AbstractPersistenceArchive_i& utStoreMapWithObjectsOfVectorWithObjects(CON_AbstractPersistenceArchive_i& ar, const std::map< Key, std::vector<T> >& m)
{
#ifdef WIN32
    ar << m.size();

    for (typename std::map< Key, std::vector<T> >::const_iterator theIterator = m.begin(); theIterator != m.end(); ++theIterator)
    {
        const typename std::map< Key, std::vector<T> >::key_type& currElement = theIterator->first;
        ar << &currElement;

        utStoreVectorOfObjectReferences<T>( ar, theIterator->second );
    }
#endif //WIN32
    return ar;
}

template< class Key, class T >
CON_AbstractPersistenceArchive_i& utLoadMapWithObjectsOfVectorWithObjects(CON_AbstractPersistenceArchive_i& ar, std::map< Key, std::vector<T> >& m)
{
#ifdef WIN32
    // make sure to pass an empty map
    ASSERT( m.empty() );

    typename std::map< Key, std::vector<T> >::size_type nSize = 0;
    ar >> nSize;

    for (typename std::map< Key, std::vector<T> >::size_type i = 0; i < nSize; ++i)
    {
        typename std::map< Key, std::vector<T> >::key_type* pKey = NULL;
#ifdef WIN32
        std::map< Key, std::vector<T> >::referent_type referent;
#else
        typename std::map< Key, std::vector<T> >::reference referent;
#endif
        ar >> pKey;
        ASSERT(pKey != NULL);

        utLoadVectorOfObjectReferences<T>( ar, referent );

        if (pKey == NULL)
            continue;

        ASSERT( m.find( *pKey ) == m.end() );
#ifdef WIN32
        m.insert( std::pair<std::map< Key, std::vector<T> >::key_type, std::map< Key, std::vector<T> >::referent_type>(*pKey, referent));
#else
           m.insert( typename std::pair<typename std::map< Key, std::vector<T> >::key_type, typename std::map< Key, std::vector<T> >::reference>(*pKey, referent));
#endif

        // insert makes a copy of pKey
        // delete object

        delete pKey;
    }
#endif //WIN32
    return ar;
}


template< class Key, class T >
CON_AbstractPersistenceArchive_i& utStoreMapOfObjectPointers(CON_AbstractPersistenceArchive_i& ar, const std::map< Key, T >& m)
{
#ifdef WIN32
    ar << m.size();

    for (typename std::map< Key, T >::const_iterator theIterator = m.begin(); theIterator != m.end(); theIterator++)
    {
        ar << theIterator->first;
        ar << theIterator->second;
    }
#endif //WIN32
    return ar;
}

template< class Key, class T >
CON_AbstractPersistenceArchive_i& utLoadMapOfObjectPointers(CON_AbstractPersistenceArchive_i& ar, std::map< Key, T >& m)
{
#ifdef WIN32
    // make sure to pass an empty map
    ASSERT( m.empty() );

    typename std::map< Key, T >::size_type nSize = 0;
    ar >> nSize;

    for (typename std::map< Key, T >::size_type i = 0; i < nSize; i++)
    {
        typename std::map< Key, T >::key_type key;
#ifdef WIN32
        std::map< Key, T >::referent_type referent = NULL;
#else
           typename std::map< Key, T >::reference referent = NULL;
#endif

        ar >> key;
        ar >> referent;

        ASSERT( m.find( key ) == m.end() );
#ifdef WIN32
        m.insert( std::pair<std::map< Key, T >::key_type, std::map< Key, T >::referent_type>(key, referent) );
#else
           m.insert( typename std::pair<typename std::map< Key, T >::key_type, typename std::map< Key, T >::reference>(key, referent) );
#endif


    }
#endif //WIN32
    return ar;
}


template< class T >
CON_AbstractPersistenceArchive_i& utStoreSet(CON_AbstractPersistenceArchive_i& ar, const std::set< T >& s)
{
#ifdef WIN32
    ar << s.size();

    for (typename std::set< T >::const_iterator theIterator = s.begin(); theIterator != s.end(); theIterator++)
    {
        ar << *theIterator;
    }
#endif //WIN32
    return ar;
}

template< class T >
CON_AbstractPersistenceArchive_i& utLoadSet(CON_AbstractPersistenceArchive_i& ar, std::set< T >& s)
{
#ifdef WIN32
    // make sure to pass an empty set
    ASSERT( s.empty() );


    typename std::set< T >::size_type nSize = 0;
    ar >> nSize;

    for (typename std::set< T >::size_type i = 0; i < nSize; i++)
    {
        typename std::set< T >::value_type CurrElement;
        ar >> CurrElement;

        s.insert(CurrElement);
    }
#endif //WIN32
    return ar;
}

//////////////////////////////////////////////////////////////////////
// CON_AbstractPersistenceArchive_i output helpers
#ifdef WIN32
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const std::pair<std::string, std::string>& Pair);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, std::pair <std::string, std::string>& Pair);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const std::vector< std::string >& String);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, std::vector< std::string >& String);


UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, SIZE size);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, SIZE& size);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, POINT point);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, POINT& point);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, RECT& rect);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, RECT& rect);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, GUID& guid);
UTILITIES_EXPORT CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, GUID& guid);
#endif

#endif // !defined(AFX_UTABSTRACTPERSISTENCEARCHIVE_I_H__42067FC7_2F96_11D3_96D4_0080C8FD5716__INCLUDED_)
