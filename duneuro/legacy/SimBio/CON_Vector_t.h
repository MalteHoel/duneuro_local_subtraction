//$29 6 july 2007 Steinstraeter  Element wise multiplication % deleted. Comments and i/o methods added.
//$28    26.03.2007  Thomas D.    Removed variable "static T errorVar" in bracket overload function
//$27   11.08.2003  TRK            added elementwise operators
//$26    11.08.2003  TRK         added mean function
//$25   31.07.2003    TRK            added exceptions
//$24    01.03.2002    Andrei        Added remove and blockremove methods
//$23    12.02.2002        TRK        blockcopy with memcpy
//$22    12.02.2002        TRK     copy function makes hardcopy always, if reset_ownmem=false
//$21    12.02.2002        TRK        reset_ownmem = false in assign function
//$20   17.01.2001  AA          #include <stdexcept> std::runtime_error and std::bad_alloc on SGI
//$19    15.11.2000  TRK            Changed Comment
//$18    14.08.2000    Frank N.    Changed terurn value of getdata() to const T*
//$17   10.07.2000  Matthias D. Removed capital letters from file names
//$16   16.06.2000  Matthias D. Made some changes to work on Linux
//$16    11.04.2000    Frank Z.    Added getMini and getMaxi()
//$15    31.03.2000    Frank N.    In CON_Vector_t(int,T=0) added default for argument
//$14    25.03.2000    TRK            Added scale function
//$13    14.03.2000    Frank N.    Handle the NOT "ownmem" case in the Serialize() function
//$12    13.03.2000    Frank N.    Worked on performance improvements
//$11    24.02.2000    Frank N.    Changed define for bad_allcoc
//$10    17.02.2000    Frank N.    Derived from CON_PersistentClass_t
//$9    20.10.1999    Frank Z.    change: CON_Vector_t :: init() reset_ownmem = true;
//$9    20.09.1999    Frank Z.    bug fix in concatenation operator
//$8    26.08.1999    Frank Z.    added sum()
//$7    20.08.1999    Frank N.    Fixed init. problem in normalise()
//$6    18.08.1999    Frank Z.    added implementation of normalise()
//$5    09.08.1999    TRK            Changed to #include <stdexcept>
//$4    01.07.1999    Frank N.    Serialization is now based on CON_AbstractPersistenceArchive_i
//$3    13.05.1999    Frank N.    changed argument name from vector to Vector in operators
//$2    10.05.1999    Frank N.    added serialization - for now for CArchive only
//$1    05.05.1999    Frank N.    split from ut_array_t.h
//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche,
//    Dr.Thomas Dierkes, Dr. Olaf Steinstraeter.
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
// File name : CON_Vector_t.h
//-------------------------------------------------------------------------

#ifndef __CON_VECTOR_T_H__
#define __CON_VECTOR_T_H__

#include <stdexcept>
#include <math.h>
#include <stdlib.h>
#include "CON_UtilitiesDef.h"
#include "Error.h"
#include "CON_AbstractPersistenceArchive_i.h"
#include "CON_PersistentClass_t.h"
// added by Steinstraeter: BGN
#include "UTI_ByteOrder_c.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
// added by Steinstraeter: END

/////////////////////////////////////////////////////////////////////////////
// some auxiliary functions

// safe square root

template<class T>
inline T root(T arg)
{
  if (arg>0)
  {
    using namespace std;
    return sqrt(arg);
  }
  else
    return 0;
}

// One-dimensional array (CON_Vector_t)
// This object can be used like an ordinary array, it is however saveguarded
// against out-of-range errors. Additional functions of e.g. CON_Vector_t algebra are
// available
// comment added by Steinstraeter: For an CON_Vector_t<T>, the array of T values are initialised
//                                 by memset(...,0, ...). Therefore, it must be possible to
//                                 initialize instances of T by sizeof(T) bytes with value 0.
//                                 Class instances are normally destroyed by this mechanism.
//------------------------------------------------------------------------------
// The following functions are available. X,Y,Z are vectors, i,j,k integers,
// a,b,c double/floats.
//
// description                        example                    comment
// ------------                       ---------------------      ---------------
// creation                           CON_Vector_t<float> X;
// creation as copy                   CON_Vector_t<double> X(Y);
// creation with allocation           CON_Vector_t<int>   X(220);
// allocation                         X.allocate(32);
//    comment added by Steinstraeter: X.allocate(0) is the same as X.clear(), this is different from CON_Matrix_t::allocate()
// assignment of data                 float a[20];X.assign(20,a);
// deallocation                       clear(X);
// test of existence                  if (exist(X))
// length of vector                   i = length(X)+length(Y);
// reference to element               a = X[23];
// assignment                         CON_Vector_t<char> X = Y;
// addition (element-wise)            X = Y + Z; Z += Y;         // operands must have equal length
// subtraction (element-wise)         X = Y - Z; Z -= Y;         // operands must have equal length
// scalar product                     double a = X * Y;          // operands must have equal length
// multiplication with scalar         X *= a; Y = Z * b;         // NOT Y = b * Z !!!
// equality / unequality              if (X == Y && Y != Z)
// Euclidian norm (yields double)     double a = norm(X);
//    comment added by Steinstraeter:   WARNING: This actually calculates the norm of vector X, but CON_Matrix_t<T>::norm() calculates
//                                      the square of the Frobenius norm !!!
// cross (vector) product             X = Y / Z;                 // operands must be of length 3
// concatination                      X = Y | Z;
// scale with another vector          X.scale(Y);                // must have same dimensions
//


// added by Steinstraeter: BGN
// ========================================================================================================================
// Some Extensions
// ===============
//
// Legend:
// -------
//
//    CON_Vector_t<T> v, u;
//    ostream out;
//
// Class information:
// ------------------
//
//    char info = v.getElementType();
//      Information about the template parameter T used for v.
//      Class method.
//      info = 'd' if T = double, info = 'f' if T = float, info = '?' otherwise.
//
// Data Access:
// ------------
//
//    CON_Vector_t<T> vec = v(i);
//      =  v[i], but without range check as long as RANGECHECKDEBUG is not set
//
// Resize:
// -------
//    v.allocate_or_clear(new_length)
//      calls clear() if new_length <= 0
//      and allocate(new_length) otherwise.
//      This is nearly the same as allocate()
//      (allocate(new_length) throws an exception
//      for new_length < 0). The major advantage of
//      this method: CON_Vector_t<T>::allocate_or_clear()
//      and CON_Matrix_t<T>::allocate_or_clear() behaves
//      equally.
//
// Data Exchange:
// --------------
//
//    v.swap(u)
//    Exchange content of v and u.
//
// I/O methods:
// ------------
//
//    out << v;
//      writes vector v to out (class ostream, e.g. cout or cerr), matlab notation
//
//    std::string s = v.toString(maxLen, setw_value, format, setprecision_value);
//    std::string s = v.toString(maxLen, setw_value, format);
//    std::string s = v.toString(maxLen, setw_value);
//    std::string s = v.toString(maxLen);
//    std::string s = v.toString();
//      converts v into string, matlab notation.
//         maxLen :              output restricted to maxLen elements,
//                               maxLen < 0: no restriction,
//                               default: -1
//         setw_value:           setw(setw_value) is used for formatting
//                               default: 0
//         format:               format of the floating point numbers:
//                                  'f' :                      setf(ios_base::fixed, ios_base::floatfield),
//                                  'e' :                      setf(ios_base::scientific, ios_base::floatfield),
//                                  otherwise, especially 'g': setf(0, ios_base::floatfield).
//                               default: 'g'.
//         setprecision_value:   setprecision(setprecision_value) is used for formating,
//                               default: -1
//
//    bool success = v.save(filename);
//      writes v to binary file filename.
//      Format for CON_Vector_t<T> and CON_Matrix_t<T>:
//         <type> <0> <0> <0> <rows> <columns> <first element> ... <last element>
//            <type> char                  : 'd' double, 'f' float
//            <0>    char                  : 0, reserved
//            <rows> int                   : number of rows
//            <columns> int                : number of columns
//            <... element> float / double : data elements row-wise
//                                           CON_Matrix_t<T> is saved as
//                                             double if T = double,
//                                             float  if T = float,
//                                             double otherwise.
//      Here:
//         <rows> = 1, <columns> = v.length(),
//         <type> = 'f' if T = float, 'd' otherwise.
//      A corresponding matlab function utMatrixLoad.m exists in <SimBioHome>/utilities/matlab.
//      success = false indicates error.
//    bool success = v.load(filename);
//      read v from binary file filename.
//      Format: see v.save(filename).
//              <rows> == 0 or <column> == 0 : v.exist() == false,
//              <rows> = 1:                   <columns> -> v.length(),
//              <columns> == 1:               <rows> -> v.length(),
//              otherwise:                    error (return value false).
//      A corresponding matlab function utMatrixSave.v exists in <SimBioHome>/utilities/matlab.
//      success = false indicates error (in this case: v is empty (v.exist() = false)).
// ========================================================================================================================
// added by Steinstraeter: END



// added by Steinstraeter: forward declaration
template<class T> class CON_Matrix_t;

template <class T>
class CON_Vector_t : public CON_PersistentClass_t< CON_Vector_t< T > >
{
  // added by Steinstraeter: allow the access of private methods by CON_Matrix_t.
  friend class CON_Matrix_t<T>;

public:
  CON_Vector_t();                                 // create zero-length CON_Vector_t
  CON_Vector_t(int,T=0);                            // create empty CON_Vector_t of certain length
  void allocate(unsigned int);                       // through away and reallocate CON_Vector_t
  void allocate_or_clear(int new_length); // added by Steinstraeter
  CON_Vector_t(const CON_Vector_t<T>&);                 // copy CON_Vector_t
  virtual ~CON_Vector_t();                                // destroy CON_Vector_t

  void clear();                             // clear CON_Vector_t
  void swap(CON_Vector_t<T> &other);    // added by Steinstraeter: exchange the content of *this and other
  void assign(int,T*);                      // assign data to CON_Vector_t
  void copy_reduced(const CON_Vector_t<T>&,const CON_Vector_t<bool>&);       // copy a selection of other CON_Vector_t
  void remove(int index);                               // remove one element
  void removeblock(int blockStart, int blockLength);       // remove a block of elements
  bool exist() const;                       // return allocation status
  unsigned length() const;                  // return size of CON_Vector_t
  T&   operator[](unsigned int);                     // reference CON_Vector_t element
  const T&   operator[](unsigned int) const;
  T&   operator()(unsigned int);        // identical to operator[] but without range check (as long as RANGECHECKDEBUG is not set)
  const T&   operator()(unsigned int) const;        // identical to operator[] but without range check (as long as RANGECHECKDEBUG is not set)
  CON_Vector_t<T>& operator=(const CON_Vector_t<T>&);         // assign CON_Vector_t to another CON_Vector_t
  CON_Vector_t<T>& operator=(const T&);                 // fill with one value
  CON_Vector_t<T>& operator+=(const CON_Vector_t<T>&);        // add another CON_Vector_t (same size)
  CON_Vector_t<T>  operator+(const CON_Vector_t<T>&);
  CON_Vector_t<T>& operator-=(const CON_Vector_t<T>&);        // subtract another CON_Vector_t (same size)
  CON_Vector_t<T>  operator-(const CON_Vector_t<T>&);
  T operator*(const CON_Vector_t<T>&) const;            // scalar multiplication with another CON_Vector_t (same size)
  CON_Vector_t<T>& operator*=(const T&);          // multiply with scalar
  CON_Vector_t<T>  operator*(const T&) const;
  bool operator==(const CON_Vector_t<T>&) const;        // equality
  bool operator!=(const CON_Vector_t<T>&) const;        // unequality
  double norm() const;                      // euclidian norm
  double mean() const;                      // average of elements
  double sum() const;                       // sum of elements
  void       normalise();                   // normalise to euclidian norm
  CON_Vector_t<T>  operator/(const CON_Vector_t<T>&) const;       // cross product (only for two 3 element vectors)
  CON_Vector_t<T>  operator|(const CON_Vector_t<T>&) const;       // concatinate both vectors
  CON_Vector_t<T>& scale(CON_Vector_t<T>&);
  // undocumented (use in emergencies only)
  T* getdata();                                                 // return pointer to data
  const T* getdata() const;                             // return pointer to data
  static char getElementType(); // added by Steinstraeter: 'f' if T = float, 'd' if T = double, '?' otherwise.
  std::string toString(int maxLen = -1,
                       int setw_value = 0,
                       char format = 'g',
                       int setprecision_value = -1) const; // added by Steinstraeter: string representation, matlab notation.
  bool save(const char *filename) const; // added by Steinstraeter: write vector elements to binary file.
  bool load(const char *filename); // added by Steinstraeter: read vector elements from binary file.

  T getMaxi() const;
  T getMini() const;

#ifdef WIN32
  // Serialization
  virtual void Serialize(class CON_AbstractPersistenceArchive_i& ar);
  friend class CON_AbstractPersistenceArchive_i& operator<<(class CON_AbstractPersistenceArchive_i& ar, const CON_Vector_t<T>& Vector);
  friend class CON_AbstractPersistenceArchive_i& operator>>(class CON_AbstractPersistenceArchive_i& ar, CON_Vector_t<T>& Vector);
#endif
private:
  void init();                         // initialise unallocated object
  void copy(const CON_Vector_t<T>&);     // common code used in assignment operator and copy constructor
  static std::ofstream *save_header(const char *filename,
                                    char outputType,
                                    int rows, int columns);  // added by Steinstraeter: open filename and save the specified header information,
                                                             //                         errors (not only i/o errors) are indicated by invalidated std::ofstream
  void rawsave(std::ofstream &f, char outputType, int start = 0, int len = -1) const;  // added by Steinstraeter: write selected vector elements to f,
                                                                                       //                         outputType = 'f' (float) for T = float,
                                                                                       //                                      'd' i(double) otherwise,
                                                                                       //                         invalid len (especially < 0) means up to the end,
                                                                                       //                         errors (not only i/o errors) are indicated by invalidated f
  static std::ifstream *read_header(const char *filename,
                                    char &fileFormat,
                                    int &rows, int &colums);  // added by Steinstraeter: open filename for reading and return header information,
                                                              //                         errors (not only i/o errors) are indicated by invalidated std::ofstream
  void rawread(std::ifstream &f,
               char fileFormat,
               int start = 0, int len = -1); // added by Steinstraeter: read len vector elements from f to *this starting with (*this)[start],
                                             //                         fileFormat must be 'f' (float) or 'd' (double),
                                             //                         invalid len (especially < 0) means up to the end (specified by length())
                                             //                         errors (not only i/o errors) are indicated by invalidated f

protected:
  // long size;                      // number of elements in CON_Vector_t
  unsigned size;                       // number of elements in CON_Vector_t
  T*       data;                       // data CON_Block_t
  bool ownmem;                         // indicates, whether the allocated memory is owned by this object and
                                       // should be copied and deleted by it
  bool reset_ownmem;                   // if set, the copy constructor resets ownmem to true
                                       // comment added by Steinstraeter:
                                       //   If src.reset_ownmem = true and src.ownmem = false, the next dest.copy(src)
                                       //   will transfer the responsibility of freeing data to dest.
                                       //   This is useful for a temporary vector that should be copied to a function caller
                                       //   by return (e.g. operator+()).

  DECLARE_PERSISTENT_CLASS( CON_Vector_t< T> )
};

// some non-element equivalents (just for more clear notation, e.g. if (norm(x)<a) instead of if (x.norm()<a)

template<class T>
inline void clear(CON_Vector_t<T>& arg)
{
  arg.clear();
}

template<class T>
inline bool exist(CON_Vector_t<T>& arg)
{
  return arg.exist();
}

template<class T> inline
unsigned length(CON_Vector_t<T>& arg)
{
  return arg.length();
}

template<class T> inline
double norm(CON_Vector_t<T>& arg)
{
  return arg.norm();
}

template<class T>
CON_Vector_t<T> :: CON_Vector_t()
{
  init();
}

// added by Steinstraeter:
//   formated stream output, matlab notation
template<class T>
std::ostream &operator<<(std::ostream &s, const CON_Vector_t<T> &v)
{
  s << "[";
  if (v.exist())
  {
    for (unsigned int i = 0; i < v.length()-1; i++)
      s << v[i] << " ";
    s << v[v.length()-1];
  }
  s << "]";

  return s;
}

// constructor and allocator: allocate CON_Vector_t with i elements of type T

template<class T>
CON_Vector_t<T> :: CON_Vector_t(int i, T a)
{
  init();
  allocate(i);
  (*this) = a;
}

template<class T>
inline void CON_Vector_t<T> :: allocate(unsigned int i)
{
  // Check for same size afterwards, if so just reuse the old memory
  if (i == size && ownmem == true && data != NULL)
  {
    ::ZeroMemory(data, sizeof(T) * size);
    return;
  }

  // Allocate data CON_Vector_t
  clear();

  if (i > 0)
  {
    try
    {
      data = new T[i];
      ownmem = true;
      ::ZeroMemory(data, sizeof(T) * i);
    }
    // in case of failure
    catch(std::bad_alloc&)
    {
      clear();
      ASSERT(FALSE);
      throw Error("CON_Vector_t::allocate(...)","Allocation failed!",ut_error_array_alloc);
    }
  }

  size = i;
}

// added by Steinstraeter: clear() if new_length <= 0, otherwise allocate().
template<class T>
inline void CON_Vector_t<T>::allocate_or_clear(int new_length)
{
  if (new_length <= 0)
    clear();
  else
    allocate(new_length);
}

// contructor: copy another CON_Vector_t

template<class T>
CON_Vector_t<T> :: CON_Vector_t(const CON_Vector_t<T>& org)
{
  init();
  copy(org);
}

// destructor and clear function: destroy CON_Vector_t and free ressources

template<class T>
CON_Vector_t<T> :: ~CON_Vector_t()
{
  clear();
}

template<class T>
inline void CON_Vector_t<T> :: clear()
{
  // in case of failure, don\B4t do anything
  if (ownmem && data)
    try
    {
      delete [] data;
    }
    catch(...)
    {
      ASSERT(FALSE);
      throw Error("CON_Vector_t::clear()","Deallocation!",ut_error_array_alloc);
    };

  init();
}

// assign array of type T to CON_Vector_t

template<class T>
inline void CON_Vector_t<T> :: assign(int i, T* orgdata)
{
  clear();
  data = orgdata;
  size = i;
  ownmem = false;
  reset_ownmem = false;       // TRK11022002, should not own the memory later either
  // this flag combination ensures "hardcopy" of object in copy operator
}

// copy selection of other CON_Vector_t

template<class T>
inline void CON_Vector_t<T> :: copy_reduced(const CON_Vector_t<T>& full,const CON_Vector_t<bool>& select)
{
  // if full CON_Vector_t is empty, just return
  if (!full.exist())
    return;

  // if selection CON_Vector_t is empty, just copy full
  if (!select.exist())
  {
    clear();
    copy(full);
    return;
  }

  // if lengths of full and selection differ
  if (select.length()!=full.length())
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::copy_reduced(...)","Wrong selection vector!",ut_error_array_selectrows);
    return;
  }

  // clear target (this) object
  clear();
  // set dimensions and memory flags
  int i,j;
  for (i=0; i<full.length(); i++)
    if (select[i]) size++;
  if (!size)
    return;

  ownmem = true;
  // try to allocate memory, if unsuccessful
  try
  {
    allocate(size);
  }
  catch(...)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::copy_reduced(...)","Allocation failed!",ut_error_array_alloc);
  };
  // copy data
  for (int i=j=0; i<full.size; i++)
    if (select[i])
      data[j] = full[j++];
}

// return allocation status

//remove
template<class T>
inline void CON_Vector_t<T> :: removeblock(int blockStart, int blockLength)
{
  //ASSERT(blockStart >= 0 && blockLength >=0 && blockLength + blockStart < size);
  if (blockStart < 0 || blockLength < 0 || blockStart+blockLength >= size)
  {
    ASSERT(false);
    throw Error("CON_Vector_t::removeblock(...)","Wrong dimensions!",ut_error_array_dimmatch);
    return;
  }

  int blockEnd = blockStart + blockLength;
  CON_Vector_t<T> target(size-blockLength);
  for (int i = 0; i < size; ++i)
    if (i < blockStart)
      target[i] = data[i];
    else
    if(i >= blockEnd)
      target[i-blockLength] = data[i];

  *this = target;
}

template<class T>
inline void CON_Vector_t<T> :: remove(int index)
{
  if (index < 0 || index >= size)
  {
    ASSERT(false);
    throw Error("CON_Vector_t::remove(...)","Wrong index!",ut_error_array_index);
    return;
  }

  CON_Vector_t<T> target(size-1);
  for (int i = 0; i < size; ++i)
    if (i != index)
      target[(i < index ? i : i - 1)] = data[i];
  //??? delete this?????
  *this = target;
}


template<class T>
inline bool CON_Vector_t<T> :: exist() const
{
  if (size && data)
    return true;
  else
    return false;
}

// return CON_Vector_t size

template<class T>
inline void CON_Vector_t<T> :: normalise()
{
  int i;
  if (!exist())
    return;
  T len = 0;
  for (i=0; i<size; i++)
    len += data[i] * data[i];
  len = sqrx(len);
  if (len == 0)
    len = 1;
  for (i=0; i<size; i++)
    data[i] /= len;
}

template<class T>
inline unsigned CON_Vector_t<T> :: length() const
{
  return size;
}

// overload brackets to reference the ith element of the CON_Vector_t

template<class T>
inline T& CON_Vector_t<T> :: operator[](unsigned int i)
{
  // static T errorVar = 0;
  if (i<0 || i>=size)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator[]","Wrong index!",ut_error_array_index);
  }
  return data[i];
}

template<class T>
inline const T& CON_Vector_t<T> :: operator[](unsigned int i) const
{
  // static T errorVar = 0;
  if (i<0 || i>=size)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator[]","Wrong index!",ut_error_array_index);
  }
  return data[i];
}

// added by Steinstraeter: exchange the content of *this and other
template<class T>
void CON_Vector_t<T>::swap(CON_Vector_t<T> &other)
{
  std::swap(size, other.size);
  std::swap(data, other.data);
  std::swap(ownmem, other.ownmem);
  std::swap(reset_ownmem, other.reset_ownmem);
}

// added by Steinstraeter: identical to operator[] but without range check as long as RANGECHECKDEBUG is not set
template<class T>
inline T& CON_Vector_t<T> :: operator()(unsigned int i)
{
#ifdef RANGECHECKDEBUG
  if (i<0 || i>=size)
  {
    std::cerr << "CON_Vector_t::operator(): " << "Wrong index!" << std::endl;
    exit(1);
  }
#endif
  return data[i];
}

// added by Steinstraeter: identical to operator[] but without range check as long as RANGECHECKDEBUG is not set
template<class T>
inline const T& CON_Vector_t<T> :: operator()(unsigned int i) const
{
#ifdef RANGECHECKDEBUG
  if (i<0 || i>=size)
  {
    std::cerr << "CON_Vector_t::operator(): " << "Wrong index!" << std::endl;
    exit(1);
  }
#endif
  return data[i];
}

// overload assignment operator

template<class T>
inline CON_Vector_t<T>& CON_Vector_t<T> :: operator=(const CON_Vector_t<T>& org)
{
  copy(org);
  return *this;
}

template<class T>
inline CON_Vector_t<T>& CON_Vector_t<T> :: operator=(const T& org)
{
  if (exist())
  {
    for (unsigned int i=0; i<size; i++)
      (* this)[i] = org;
  }
  return *this;
}

// addition

template<class T>
inline CON_Vector_t<T>& CON_Vector_t<T> :: operator+=(const CON_Vector_t<T>& org)
{
  // test size of the source CON_Vector_t
  if (org.size!=size || !org.data || !data)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator+=","Mismatch between operands!",ut_error_array_dimmatch);
    return *this;
  }
  // copy data
  for (int i=0; i<org.size; i++)
    data[i] += org[i];

  return *this;
}

template<class T>
inline CON_Vector_t<T> CON_Vector_t<T> :: operator+(const CON_Vector_t<T>& org)
{
  // create target object
  CON_Vector_t<T> target;
  // test size of the source CON_Vector_t
  if (org.size!=size || !org.data || !data)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator+","Mismatch between operands!",ut_error_array_dimmatch);
    return target;
  }
  try
  {
    target.allocate(size);
  }
  catch (...)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator+","Allocation failed!",ut_error_array_alloc);
    return CON_Vector_t<T>();
  }
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  // copy data
  for (int i=0; i<org.size; i++)
    target[i] = data[i] + org[i];

  return target;
}

// subtraction

template<class T>
inline CON_Vector_t<T>& CON_Vector_t<T> :: operator-=(const CON_Vector_t<T>& org)
{
  // test size of the source CON_Vector_t
  if (org.size!=size || !org.data || !data)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator-=","Mismatch between operands!",ut_error_array_dimmatch);
    return *this;
  }
  // copy data
  for (int i=0; i<org.size; i++)
    data[i] -= org[i];

  return *this;
}

template<class T>
inline CON_Vector_t<T> CON_Vector_t<T> :: operator-(const CON_Vector_t<T>& org)
{
  // create target object
  CON_Vector_t<T> target;
  // test size of the source CON_Vector_t
  if (org.size!=size || !org.data || !data)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator-","Mismatch between operands!",ut_error_array_dimmatch);
    return target;
  }
  try
  {
    target.allocate(size);
  }
  catch (...)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator-","Allocation failed!",ut_error_array_alloc);
    return CON_Vector_t<T>();
  }
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  // copy data
  for (int i=0; i<org.size; i++)
    target[i] = data[i] - org[i];

  return target;
}

// scalar multiplication

template<class T>
inline T CON_Vector_t<T> :: operator*(const CON_Vector_t<T>& org) const
{
  // create target object
  T target = 0;
  // test size of the source CON_Vector_t
  if (org.size!=size || !org.data || !data)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator*","Mismatch between operands!",ut_error_array_dimmatch);
    return target;
  }
  // carry out multiplication
  for (int i=0; i<org.size; i++)
    target += data[i] * org[i];

  return target;
}

// multiplication with scalar

template<class T>
inline CON_Vector_t<T>& CON_Vector_t<T> :: operator*=(const T& org)
{
  if (exist())
  {
    for (int i=0; i<size; i++)
      data[i] *= org;
  }

  return *this;
}

template<class T>
inline CON_Vector_t<T> CON_Vector_t<T> :: operator*(const T& org) const
{
  // create target object
  CON_Vector_t<T> target;
  try
  {
    target.allocate(size);
  }
  catch (...)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator*","Allocation failed!",ut_error_array_alloc);
    return CON_Vector_t<T>();
  }
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  if (exist())
  {
    for (int i=0; i<size; i++)
      target[i] = data[i] * org;
  }

  return target;
}

// equality

template<class T>
inline bool CON_Vector_t<T> :: operator==(const CON_Vector_t<T>& org) const
{
  // if both vectors are empty, the ARE equal
  if (!exist() && !org.exist())
    return true;
  // if the vectors have different lengths, the are not equal
  if (size != org.size)
    return false;
  // check for equality of the arguments
  for (int i=0; i<size; i++)
    if (data[i] != org[i])
      return false;

  return true;
}

template<class T>
inline bool CON_Vector_t<T> :: operator!=(const CON_Vector_t<T>& org) const
{
  return !((*this)==org);
}

// euclidian norm

template<class T>
inline double CON_Vector_t<T> :: norm() const
{
  // if CON_Vector_t is empty, return 0
  if (!exist())
    return 0.0;
  // multiply this CON_Vector_t with itself and return (safe) square root of it
  double N=0;
  for (int i=0; i<size; i++)
    N += (double)data[i] * (double)data[i];
  return root(N);
}

template<class T>
inline double CON_Vector_t<T> :: sum() const
{
  // if CON_Vector_t is empty, return 0
  if (!exist())
    return 0.0;
  // add elements of CON_Vector_t
  double N=0;
  for (int i=0; i<size; i++)
    N += data[i];
  return (N);
}

template<class T>
inline double CON_Vector_t<T> :: mean() const
{
  // if CON_Vector_t is empty, return 0
  if (!exist())
    return 0.0;
  return (sum()/size);
}

// cross product (only for two 3 element vectors)

template<class T>
inline CON_Vector_t<T> CON_Vector_t<T> :: operator/(const CON_Vector_t<T>& org) const
{
  // create target object
  CON_Vector_t<T> target;
  // test size of the source and this CON_Vector_t (must be 3)
  if (org.size!=3 || size!=3 || !org.data || !data)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator/","Mismatch between operands!",ut_error_array_dimmatch);
    return target;
  }
  try
  {
    target.allocate(size);
  }
  catch (...)
  {
    ASSERT(FALSE);
    throw Error("CON_Vector_t::operator/","Allocation failed!",ut_error_array_alloc);
    return CON_Vector_t<T>();
  }
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  // carry out multiplication
  target[0] = data[1]*org[2]-data[2]*org[1];
  target[1] = data[2]*org[0]-data[0]*org[2];
  target[2] = data[0]*org[1]-data[1]*org[0];

  return target;
}

// concatination of two vectors

template<class T>
inline CON_Vector_t<T> CON_Vector_t<T> :: operator|(const CON_Vector_t<T>& org) const
{
  // create target object
  CON_Vector_t<T> target;
  if (exist() || org.exist())
  {
    try
    {
      target.allocate(size+org.size);
    }
    catch (...)
    {
      ASSERT(FALSE);
      throw Error("CON_Vector_t::operator|","Allocation failed!",ut_error_array_alloc);
      return CON_Vector_t<T>();
    }
    // make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
    // fill in values
    int i,j;
    for (i=j=0; i<size; i++,j++) target[j] = data[i];
    for (i=0; i<org.size; i++,j++)
      target[j] = org[i];
  }
  return target;
}

template<class T> CON_Vector_t<T>& CON_Vector_t<T> :: scale(CON_Vector_t<T>& scalor)
{
  if (scalor.length() == length())
    for (int i=0; i<length(); i++)
      data[i] *= scalor.data[i];
  return *this;
}

// harakiri function: do not use unless authorised
template<class T>
inline T* CON_Vector_t<T> :: getdata()
{
  return data;
}

// harakiri function: do not use unless authorised
template<class T>
inline const T* CON_Vector_t<T> :: getdata() const
{
  return data;
}

// added by Steinstraeter: 'f' if T = float, 'd' if T = double, '?' otherwise.
//                         specialized version: double
template<>
inline char CON_Vector_t<double>::getElementType()
{
  return 'd';
}

// added by Steinstraeter: 'f' if T = float, 'd' if T = double, '?' otherwise.
//                         specialized version: float
template<>
inline char CON_Vector_t<float>::getElementType()
{
  return 'f';
}

// added by Steinstraeter: 'f' if T = float, 'd' if T = double, '?' otherwise.
//                         standard version
template<class T>
inline char CON_Vector_t<T>::getElementType()
{
  return '?';
}

// added by Steinstraeter:
//   string conversion, matlab notation
template<class T>
std::string CON_Vector_t<T>::toString(int maxLen, int setw_value, char format, int setprecision_value) const
{
  std::ostringstream s;

  if (format == 'f')
    s << std::fixed;
  else if (format == 'e')
    s << std::scientific;

  s << std::setprecision(setprecision_value);

  s << "[";
  if (exist())
  {
    if (maxLen < 0 || maxLen > length())
      maxLen = length();

    if (maxLen > 0)
    {
      for (int i = 0; i < maxLen-1; i++)
        s << std::setw(setw_value) << (*this)[i] << " ";
      s << std::setw(setw_value) << (*this)[maxLen-1];
    }

    if (maxLen < length())
    {
      int furtherElements = length()-maxLen;
      char pl[] = "s";
      if (furtherElements == 1)
        pl[0] = 0;
      s << " ... <" << furtherElements << " further element" << pl << ">";
    }
  }
  s << "]";

  return s.str();
}

// added by Steinstraeter: save vector elements to binary file
template<class T>
inline bool CON_Vector_t<T>::save(const char *filename) const
{
  // number of elements converted at once
  const int packageSize = 1048576; // 1MB

  char outputFormat = 'd';
  if (getElementType() == 'f')
    outputFormat = 'f';
  std::ofstream *fptr = save_header(filename, outputFormat, 1, length());

  if (*fptr)
  {
    if (getElementType() == 'f' || getElementType() == 'd')
      rawsave(*fptr, outputFormat);
    else
    {
      int pos = 0;
      while (pos < size && (*fptr))
      {
        rawsave(*fptr, outputFormat, pos, packageSize);
        pos += packageSize;
      }
    }
  }

  bool success = (*fptr?true:false);
  fptr->close();
  delete fptr;
  return success;
}

// added by Steinstraeter: read vector elements from binary file
template<class T>
inline bool CON_Vector_t<T>::load(const char *filename)
{
  // number of elements converted at once
  const int packageSize = 1048576; // 1MB

  char fileFormat;
  int rows, columns;
  std::ifstream *fptr = read_header(filename, fileFormat, rows, columns);

  int len = -1;
  if (rows * columns == 0)
    len = 0;
  else if (rows == 1 && columns > 0)
    len = columns;
  else if (rows > 0 && columns == 1)
    len = rows;

  if ((fileFormat == 'd' || fileFormat == 'f') && len >= 0)
  {
    allocate_or_clear(len);
    if (getElementType() == fileFormat)
      rawread(*fptr, fileFormat);
    else
    {
      int pos = 0;
      while (pos < size && (*fptr))
      {
        rawread(*fptr, fileFormat, pos, packageSize);
        pos += packageSize;
      }
    }
  }
  else
    fptr->setstate(std::ios_base::badbit);

  bool success = (*fptr?true:false);
  fptr->close();
  delete fptr;
  if (!success)
    clear();
  return success;
}


//------------------------------------------------------------------------------
// private functions

// initialise unallocated object

template<class T>
inline void CON_Vector_t<T> :: init()
{
  data = NULL;
  size = 0;
  ownmem = false;
  reset_ownmem = true;
}

// common code used in assignment operator and copy constructor

template<class T>
inline void CON_Vector_t<T> :: copy(const CON_Vector_t<T>& org)
{
  // clear target (this) object
  clear();
  // copy dimensions and memory flags
  size = org.size;
  if (org.reset_ownmem)
    ownmem = true;
  else
    ownmem = org.ownmem;
  if (org.ownmem || !org.reset_ownmem)       // TRK11022002, if not ownmem=false AND reset_ownmem=true,
  // as used in temporary objects that are put on stack by operator functions
  {
    // try to allocate memory, if unsuccessful
    try
    {
      allocate(org.size);
    }
    catch (...)
    {
      ASSERT(FALSE);
      throw Error("CON_Vector_t::copy()","Allocation failed!",ut_error_array_alloc);
      return;
    }
    // copy data
    memcpy(data,org.data,size*sizeof(T));              //TRK11022002
    ownmem = true;                                                         //TRK11022002

    // for (int i=0;i<org.size;i++)
    //  data[i] = org[i];
  }
  else       // if the data is not owned by the original, we do not copy them, but just reference
  {
    data = org.data;
  }
}

template <class T>
T CON_Vector_t<T>::getMaxi() const
{
  T tMaxi = 0;
  if (exist())
  {
    tMaxi = data[0];
    for (long nSample = 1; nSample < size; ++nSample)
      tMaxi = max(tMaxi, data[nSample]);
  }

  return tMaxi;
}

template <class T>
T CON_Vector_t<T>::getMini() const
{
  T tMini = 0;
  if (exist())
  {
    tMini = data[0];
    for (long nSample = 1; nSample < size; ++nSample)
      tMini = min(tMini, data[nSample]);
  }

  return tMini;
}

// added by Steinstraeter: open filename and save the specified header information
template<class T>
inline std::ofstream *CON_Vector_t<T>::save_header(const char *filename, char outputType, int rows, int columns)
{
  std::ofstream *fptr = new std::ofstream;

  if (outputType != 'd' && outputType != 'f')
  {
    fptr->setstate(std::ios_base::badbit);
    return fptr;
  }

  if (rows < 0)
    rows = 0;

  if (columns < 0)
    columns = 0;

  if (filename)
  {
    fptr->open(filename, std::ios_base::out | std::ios_base::binary |std::ios_base::trunc);
    if (*fptr)
    {
      char code[4] = {outputType, UTI_ByteOrder_c::get_host_byte_order(), sizeof(int), static_cast<char>(outputType == 'd' ? sizeof(double) : sizeof(float))};
      fptr->write(reinterpret_cast<const char*>(code), sizeof(code));
      if (!(*fptr))
        return fptr;
      int dims[] = {rows, columns};
      fptr->write(reinterpret_cast<const char*>(dims), sizeof(dims));
    }
  }
  else
    fptr->setstate(std::ios_base::badbit);

  return fptr;
}

// added by Steinstraeter: write the selected vector elements to f
template <class T>
inline void CON_Vector_t<T>::rawsave(std::ofstream &f, char outputType, int start, int len) const
{
  if (!f)
    return;

  if (start >= size)
    return;

  int maxlen = size - start;
  if (len < 0 || len > maxlen)
    len = maxlen;

  if (len)
  {
    if (getElementType() == outputType)
      f.write(reinterpret_cast<const char *>(data + start), len*sizeof(T));
    else
    if (outputType == 'd')
    {
      double *dataConverted = new double[len];
      for (int i = 0; i < len; i++)
        dataConverted[i] = static_cast<double>(data[i+start]);
      f.write(reinterpret_cast<const char *>(dataConverted), len*sizeof(double));
      delete [] dataConverted;
    }
    else
      f.setstate(std::ios_base::badbit);
  }
}

// added by Steinstraeter: open filename for reading and return header information
template<class T>
std::ifstream *CON_Vector_t<T>::read_header(const char *filename, char &fileFormat, int &rows, int &colums)
{
  std::ifstream *fptr = new std::ifstream;

  if (filename)
  {
    fptr->open(filename, std::ios_base::in | std::ios_base::binary);
    if (*fptr)
    {
      char code[4];
      fptr->read(code, sizeof(code));
      fileFormat = code[0];
      char byteOrder = UTI_ByteOrder_c::get_host_byte_order();
      if (   (fileFormat != 'f' && fileFormat != 'd')
             || (code[1] != '?' && byteOrder != '?' && code[1] != byteOrder)
             || code[2] != sizeof(int)
             || code[3] != (fileFormat == 'd' ? sizeof(double) : sizeof(float)))
        fptr->setstate(std::ios_base::badbit);

      if (!(*fptr))
        return fptr;
      int dims[2];
      fptr->read(reinterpret_cast<char *>(dims), sizeof(dims));
      rows = dims[0];
      colums = dims[1];
    }
  }
  else
    fptr->setstate(std::ios_base::badbit);

  return fptr;
}

// added by Steinstraeter: read vector elements from f
template<class T>
void CON_Vector_t<T>::rawread(std::ifstream &f, char fileFormat, int start, int len)
{
  if (!f)
    return;

  if (start >= size)
    return;

  int maxlen = size - start;
  if (len < 0 || len > maxlen)
    len = maxlen;

  if (len)
  {
    if (getElementType() == fileFormat)
      f.read(reinterpret_cast<char *>(data + start), len*sizeof(T));
    else
    if (fileFormat == 'd')
    {
      double *rawdata = new double[len];
      f.read(reinterpret_cast<char *>(rawdata), len*sizeof(double));
      for (int i = 0; i < len; i++)
        data[i+start] = static_cast<T>(rawdata[i]);
      delete [] rawdata;
    }
    else
      f.setstate(std::ios_base::badbit);
  }
}


/////////////////////////////////////////////////////////////////////////////
// CON_Vector_t Serialization
//
#ifdef WIN32
template<class T>
inline void CON_Vector_t<T>::Serialize(class CON_AbstractPersistenceArchive_i& ar)
{
  if (ar.IsStoring())
  {
    // Check for the "ownmem" flag:
    if (exist() && !ownmem)
    {
      ASSERT(false);
#ifdef _AFXDLL
      AfxThrowArchiveException(CArchiveException::badClass, _T(""));
#else
      throw std::runtime_error(_T("badClass"));
#endif
    }

    ar << size;
    ar.Write(data, size * sizeof(T));
  }
  else
  {
    ar >> size;

    allocate(size);
    reset_ownmem = false;

    ar.Read(data, size * sizeof(T));
  }
}

template<class T>
inline CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const CON_Vector_t<T>& Vector)
{
  const_cast< CON_Vector_t<T>& >(Vector).Serialize(ar);

  return ar;
}

template<class T>
inline CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, CON_Vector_t<T>& Vector)
{
  Vector.Serialize(ar);

  return ar;
}
#endif
#endif // __utVector_t_H__
