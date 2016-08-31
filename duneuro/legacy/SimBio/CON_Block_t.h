//$20   20.07.2003 TRK          added getData() (Harakiri-Funktion)
//$19    21.03.2002 Anwander A.  size was not set in the function assign
//$18    12.03.2002 Anwander A.  changed for (j=k=0;j<layers;j++,k++) to for (j=0,k=0;j<layers;j++,k++)
//$17    12.02.2002        TRK        blockcopy with memcpy
//$16    12.02.2002        TRK     copy function makes hardcopy always, if reset_ownmem=false
//$15    12.02.2002        TRK        reset_ownmem = false in assign function
//$14   02.07.2002  Matthias D. Adapted for Linux
//$13    10.05.2001  Maurice        Changed operator unsigned short() to operator unsigned short()     const
//$12    08.05.2001  Maurice        Added casting operator: operator unsigned short();
//$11    18.01.2001    Frank N.    Added getMini and getMaxi()
//$10   10.07.2000     Matthias D.    Removed Capital letters from file names
//$9    16.06.2000     Matthias D     Made some changes to work on Linux
//$8    14.03.2000    Frank N.    1. Worked on performance improvements 2. Handle the NOT "ownmem" case in the Serialize() function
//$7    24.02.2000    Frank N.    Changed define for bad_allcoc
//$6    21.02.2000    Frank Z.    bug fix in allocation routine
//$5    17.02.2000    Frank N.    Derived from CON_PersistentClass_t
//$4    09.08.1999    TRK            Changed to #include <stdexcept>
//$3    01.07.1999    Frank N.    Serialization is now based on CON_AbstractPersistenceArchive_i
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
// File name : CON_Block_t.h
//-------------------------------------------------------------------------

#ifndef __CON_BLOCK_T_H__
#define __CON_BLOCK_T_H__

#include "CON_UtilitiesDef.h"

#include "Error.h"

#include "CON_AbstractPersistenceArchive_i.h"
#include "CON_PersistentClass_t.h"

// Three-dimensional array (CON_Block_t)
// This object can be used like and ordinary array, it is however saveguarded
// against out-of-range errors.
//------------------------------------------------------------------------------
// The following functions are available. X,Y,Z are blocks, i,j,k,rows,cols
// integers, a,b,c double/floats.
//
// description                        example                    comment
// ------------                       ---------------------      ---------------
// creation                           CON_Block_t<float> X;
// creation as copy                   CON_Block_t<long> X(Y);
// creation and allocation            CON_Block_t<char> X(lays,rows,cols);
// allocation                         X.allocate(lays,rows,cols)
// assignment of data                 int i[40];X.assign(2,2,10,i);
// deallocation                       clear(X);
// test of existence                  if (exist(X))             // true if CON_Matrix_t is allocated with both dimensions > 0
// dimensions of CON_Matrix_t               i = length(X)*height(X);  // length is number of columns, height number of rows
//                                    depth(X);                 // depth is the number of layers
// ref. to layers, rows, elements     CON_Matrix_t<int> Z = X[3];     // reference to e.g. column vectors not possible
//                                    CON_Vector_t<int> z = X[3][2];
//                                    int a = X[3][1][4];
// ref. to alternative layers         X.layer_constheight(3);
//                                    X.layer_constlength(22);
// overwrite layers of diff. orient.  X.setlayer_constdepth(X,22); // dimensions must match
//                                    X.setlayer_constheight(Y,2);
//                                    X.setlayer_constlength(Z,0);
// assignment                         CON_Matrix_t<char> X = Y;
// addition (element-wise)            X = Y + Z; Z += Y;         // operands must have equal dimensions
// subtraction (element-wise)         X = Y - Z; Z -= Y;         // operands must have equal dimensions
// multiplication with scalar         X *= a; Y = Z * b;         // NOT Y = b * Z !!!
// equality / unequality              if (X == Y && Y != Z)
// transpose                          transpose_rowcol(X);       // NO assignment (instead of Y=transpose(X), use Y=X;transpose(Y)),
//                                    transpose_rowlay(X);       // CAUTION: CON_Matrix_t is manipulated
//                                    transpose_collay(X);       // CAUTION: in asymmetric case the CON_Matrix_t is temporarily copied
//                                                               //
// mirror                             mirror_horizontal(X);      // reversing column indices
//                                    mirror_vertical(Z);        // reversing row indices
//                                    mirror_depth(Y);           // reversing layer indices
//                                                               // CAUTION: as with transpose, no assignment is possible. This
//                                                               //          is a safety measure, because the CON_Matrix_t is manipulated
// concatination                      X = ((Y | Z) & W ^ V)      // first concatinates columns of Y and Z, then concatinates the rows of the result
//                                                               // with the rows of W, and finally concatinate the layers to the layers of V,
//                                                               // dimensions must match.
// casting unsigned short              operator unsigned short(); // cast from any type to unsigned short
//------------------------------------------------------------------------------

template <class T>
class CON_Block_t : public CON_PersistentClass_t< CON_Block_t< T > >
{
public:
    CON_Block_t();                            // create CON_Block_t with 0x0x0 elements
    CON_Block_t(int,int,int,T=0);             // create empty CON_Block_t of given dimensions
    void allocate(int,int,int);         // through away and reallocate CON_Block_t
    CON_Block_t(const CON_Block_t<T>&);             // copy CON_Block_t
    virtual ~CON_Block_t();                           // delete CON_Block_t

    void clear();                       // clear CON_Block_t
    void assign(int,int,int,T*);        // assign data to CON_Block_t
    bool exist() const;                 // return allocation status
    unsigned length() const;            // horizontal dimension (columns)
    unsigned height() const;            // vertical dimension   (rows)
    unsigned depth()  const;            // third dimension      (layers)
    CON_Matrix_t<T>& operator[](int);         // reference a CON_Block_t layer, can be further de-referenced using the CON_Matrix_t operator
    CON_Matrix_t<T>& operator[](int) const;
    CON_Matrix_t<T>& layer_constheight(int);  // access alternative layers (involves copying)
    CON_Matrix_t<T>& layer_constlength(int);
    void       setlayer_constdepth(const CON_Matrix_t<T>&, int);
    void       setlayer_constheight(const CON_Matrix_t<T>&, int);
    void       setlayer_constlength(const CON_Matrix_t<T>&, int);
    CON_Block_t<T>& operator=(const CON_Block_t<T>&);   // assign CON_Block_t to another CON_Block_t
    CON_Block_t<T>& operator=(const T&);          // fill with one value
    CON_Block_t<T>& operator+=(const CON_Block_t<T>&);  // add another CON_Block_t (same size)
    CON_Block_t<T>  operator+(const CON_Block_t<T>&);
    CON_Block_t<T>& operator-=(const CON_Block_t<T>&);  // subtract another CON_Block_t (same size)
    CON_Block_t<T>  operator-(const CON_Block_t<T>&);
    CON_Block_t<T>& operator*=(const T&);     // multiply with scalar
    CON_Block_t<T>  operator*(const T&)    const;
    bool operator==(const CON_Block_t<T>&) const; // equality
    bool operator!=(const CON_Block_t<T>&) const; // unequality
    void transpose_rowcol();            // swap rows and columns
    void transpose_rowlay();            // swap rows and layers
    void transpose_collay();            // swap colums and layers
    void mirror_horizontal();           // X[i][j][k] becomes X[i][j][p-k]
    void mirror_vertical();             // X[i][j][k] becomes X[i][n-j][k]
    void mirror_depth();                // X[i][j][k] becomes X[m-i][j][k]
    CON_Block_t<T>  operator|(const CON_Block_t<T>&);   // concatinate horizontally
    CON_Block_t<T>  operator&(const CON_Block_t<T>&);   // concatinate vertically
    CON_Block_t<T>  operator^(const CON_Block_t<T>&);   // concatinate depth
    T* getdata();                                            // return pointer to data
    const T* getdata() const;                                // return pointer to data


// casting operators
    operator unsigned short() const;

    T getMaxi() const;
    T getMini() const;

// Serialization
#ifdef WIN32
    virtual void Serialize(CON_AbstractPersistenceArchive_i& ar);
    friend CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const CON_Block_t<T>& block);
    friend CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, CON_Block_t<T>& block);
#endif
private:
   void init();                                    // initialise unallocated object
   void copy(const CON_Block_t<T>&);                // common code used in assignment operator and copy constructor

protected:
   unsigned                layers,rows,columns;    // dimensions of CON_Block_t
   unsigned                size;                   // size of data allocated CON_Block_t (0 if no allocation)
   CON_Matrix_t<T>*        datalayers;             // CON_Matrix_t objects for layers
   T*                    data;                   // data CON_Block_t
   bool                    ownmem;                 // indicates, whether the allocated memory is owned by this object and
                                                // should be copied and deleted by it
   bool                    reset_ownmem;           // if set, the copy constructor resets ownmem to true
   CON_Vector_t<unsigned> hashtable;                // indicates offsets of the layer pointers

   DECLARE_PERSISTENT_CLASS( CON_Block_t< T> )
};

// some non-element equivalents (just for more clear notation, e.g. if (exist(x)) instead of if (x.exist())

template<class T> inline void     clear (CON_Block_t<T>& arg) {arg.clear();}
template<class T> inline bool     exist (CON_Block_t<T>& arg) {return arg.exist();}
template<class T> inline unsigned length(CON_Block_t<T>& arg) {return arg.length();}
template<class T> inline unsigned height(CON_Block_t<T>& arg) {return arg.height();}
template<class T> inline unsigned depth (CON_Block_t<T>& arg) {return arg.depth();}
template<class T> inline void     transpose_rowcol (CON_Block_t<T>& arg) {arg.transpose_rowcol();}
template<class T> inline void     transpose_rowlay (CON_Block_t<T>& arg) {arg.transpose_rowlay();}
template<class T> inline void     transpose_collay (CON_Block_t<T>& arg) {arg.transpose_collay();}
template<class T> inline void     mirror_horizontal(CON_Block_t<T>& arg) {arg.mirror_horizontal();}
template<class T> inline void     mirror_vertical  (CON_Block_t<T>& arg) {arg.mirror_vertical();}
template<class T> inline void     mirror_depth     (CON_Block_t<T>& arg) {arg.mirror_depth();}

// constructor: just initialise

template<class T>
CON_Block_t<T> :: CON_Block_t()
 {
  init();
 }

// constructor and allocator: allocate CON_Block_t with i x j x k elements of type T

template<class T>
CON_Block_t<T> :: CON_Block_t(int i,int j,int k,T a)
 {
  init();
  allocate(i,j,k);
  (*this) = a;
 }

template<class T>
inline void CON_Block_t<T> :: allocate(int i,int j,int k)
 {
  if (i<0 || j<0 || k<0) throw Error("CON_Block_t::allocate(...)","Only positive dimensions allowed.",ut_error_array_index);
  if (i==0 || j==0 || k==0) return;
  clear();
  // allocate data CON_Block_t
  try
   {
    data = new T[i*j*k];
    datalayers = new CON_Matrix_t<T>[i];
    hashtable.allocate(i);
    for (int ii=0;ii<i;ii++) {datalayers[ii].assign(j,k,data+ii*j*k);hashtable[ii]=ii*j*k;}
    ownmem = true;
   }
  // in case of failure, throw exception and leave
  catch (std::bad_alloc) // seems to be Borland specific (otherwise std::bad_alloc)
   {
    clear();
    throw Error("CON_Block_t::allocate(...)","Memory could not be allocated.",ut_error_array_alloc);
   }
  layers = i;
  rows = j;
  columns = k;
  size = i*j*k;
 }

// contructor: copy another CON_Block_t

template<class T>
CON_Block_t<T> :: CON_Block_t(const CON_Block_t<T>& org)
 {
  init();
  copy(org);
 }

// destructor and clear function: destroy CON_Block_t and free ressources

template<class T>
CON_Block_t<T> :: ~CON_Block_t()
 {
  clear();
 }

template<class T>
inline void CON_Block_t<T> :: clear()
 {
  // in case of failure, dont do anything
  if (datalayers)     {try {delete [] datalayers;} catch(...) {};}
  if (data && ownmem) {try {delete [] data;} catch(...) {};}
  init();
 }

// assign array of type T to CON_Block_t

template<class T>
inline void CON_Block_t<T> :: assign(int i,int j,int k, T* orgdata)
 {
  if (i<0 || j<0 || k<0) throw Error("CON_Block_t::assign(...)","Only positive dimensions allowed.",ut_error_array_index);
  if (i==0 || j==0 || k==0) return;
  clear();
  data = orgdata;
  try
   {
    datalayers = new CON_Matrix_t<T>[i];
    hashtable.allocate(i);
    for (int ii=0;ii<i;ii++) {datalayers[ii].assign(j,k,data+ii*j*k);hashtable[ii]=ii*j*k;}
   }
  catch (std::bad_alloc) // seems to be Borland specific (otherwise std::bad_alloc)
   {
    clear();
    throw Error("CON_Block_t::assign(...)","Memory could not be allocated.",ut_error_array_alloc);
   }
  layers = i;
  rows = j;
  columns = k;
  size = i * j *k;          // Anwander A. 21/03/2002 size was not set
  ownmem = false;
  reset_ownmem = false; // TRK11022002, should not own the memory later either
 }

// return allocation status

template<class T>
inline bool CON_Block_t<T> :: exist() const
 {
  if (rows && columns && layers && data && datalayers) return true;
  else                                     return false;
 }

// return dimensions

template<class T>
inline unsigned CON_Block_t<T> :: length() const
 {
  return columns;
 }

template<class T>
inline unsigned CON_Block_t<T> :: height() const
 {
  return rows;
 }

template<class T>
inline unsigned CON_Block_t<T> :: depth() const
 {
  return layers;
 }

// overload brackets to reference the ith layer of the CON_Block_t

template<class T>
inline CON_Matrix_t<T>& CON_Block_t<T> :: operator[](int i)
 {
  if (i<0 || i>=layers)
   throw Error("CON_Block_t::operator[]","Index out of range.",ut_error_array_range);
  return datalayers[i];
 }

template<class T>
inline CON_Matrix_t<T>& CON_Block_t<T> :: operator[](int i) const
 {
  if (i<0 || i>=layers)
   throw Error("CON_Block_t::operator[]","Index out of range.",ut_error_array_range);
  return datalayers[i];
 }

// overload assignment operator

template<class T>
inline CON_Block_t<T>& CON_Block_t<T> :: operator=(const CON_Block_t<T>& org)
 {
  copy(org);
  return *this;
 }

template<class T>
inline CON_Block_t<T>& CON_Block_t<T> :: operator=(const T& org)
 {
  if (exist())
   {
    int i,j,k;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (* this)[i][j][k] = org;
   }
  return *this;
 }

// define addition and substraction operators

template<class T>
inline CON_Block_t<T>& CON_Block_t<T> :: operator+=(const CON_Block_t<T>& org)
 {
  // test size of the source CON_Block_t
  if (org.rows!=rows || org.columns!=columns || org.layers!=layers || !org.data || !data)
      throw Error("CON_Block_t::operator+=","Dimensions of operands differ.",ut_error_array_dimmatch);
  // copy data
  try
   {
    int i,j,k;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (*this)[i][j][k] += org[i][j][k];
   }
  catch (...) {throw;}
  return *this;
 }

template<class T>
inline CON_Block_t<T> CON_Block_t<T> :: operator+(const CON_Block_t<T>& org)
 {
  // test size of the source CON_Block_t
  if (org.rows!=rows || org.columns!=columns || org.layers!=layers || !org.data || !data)
      throw Error("CON_Block_t::operator+","Dimensions of operands differ.",ut_error_array_dimmatch);
  // create target object
  CON_Block_t<T> target;
  try {target.allocate(layers,rows,columns);}
  catch (...) {throw;}
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  // copy data
  try
   {
    int i,j,k;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) target[i][j][k] = (*this)[i][j][k] + org[i][j][k];
   }
  catch (...) {throw;}
  return target;
 }

template<class T>
inline CON_Block_t<T>& CON_Block_t<T> :: operator-=(const CON_Block_t<T>& org)
 {
  // test size of the source CON_Block_t
  if (org.rows!=rows || org.columns!=columns || org.layers!=layers || !org.data || !data)
      throw Error("CON_Block_t::operator-=","Dimensions of operands differ.",ut_error_array_dimmatch);
  // copy data
  try
   {
    int i,j,k;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (*this)[i][j][k] -= org[i][j][k];
   }
  catch (...) {throw;}
  return *this;
 }

template<class T>
inline CON_Block_t<T> CON_Block_t<T> :: operator-(const CON_Block_t<T>& org)
 {
  // test size of the source CON_Block_t
  if (org.rows!=rows || org.columns!=columns || org.layers!=layers || !org.data || !data)
      throw Error("CON_Block_t::operator-","Dimensions of operands differ.",ut_error_array_dimmatch);
  // create target object
  CON_Block_t<T> target;
  try {target.allocate(layers,rows,columns);}
  catch (...) {throw;}
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  // copy data
  try
   {
    int i,j,k;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) target[i][j][k] = (*this)[i][j][k] - org[i][j][k];
   }
  catch (...) {throw;}
  return target;
 }

// multiplication with scalar

template<class T>
inline CON_Block_t<T>& CON_Block_t<T> :: operator*=(const T& org)
 {
  if (exist())
   try
    {
     int i,j,k;
     for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (*this)[i][j][k] *= org;
    }
   catch (...) {throw;}
  return *this;
 }

template<class T>
inline CON_Block_t<T> CON_Block_t<T> :: operator*(const T& org) const
 {
  // create target object
  CON_Block_t<T> target;
  try {target.allocate(layers,rows,columns);}
  catch (...) {throw;}
  // make sure, the copy constructor leaves the allocated memory alone
  target.ownmem = false;
  target.reset_ownmem = true;
  // copy data
  try
   {
    int i,j,k;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) target[i][j][k] = (*this)[i][j][k] * org;
   }
  catch (...) {throw;}
  return target;
 }

// equality

template<class T>
inline bool CON_Block_t<T> :: operator==(const CON_Block_t<T>& org) const
 {
  // if both blocks are empty, the ARE equal
  if (!exist() && !org.exist()) return true;
  // if the blocks have different dimensions, the are not equal
  if (layers != org.layers) return false;
  if (rows != org.rows) return false;
  if (columns != org.columns) return false;
  // check for equality of the elements
  int i,j,k;
  try {for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) if ((*this)[i][j][k] != org[i][j][k]) return false;}
  catch (...) {throw;}
  return true;
 }

template<class T>
inline bool CON_Block_t<T> :: operator!=(const CON_Block_t<T>& org) const
 {
  return !((*this)==org);
 }

// transpose CON_Block_t, exchange rows and columns

template<class T> void CON_Block_t<T> :: transpose_rowcol()
 {
  if (rows==columns) // in the symmetric case, we can transpose without actually copying the data
   {
    int i,j,k;
    T   a;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<j;k++) {a=(*this)[i][j][k];(*this)[i][j][k]=(*this)[i][k][j];(*this)[i][k][j]=a;}
   }
  else // for the unsymmetric case an intermediate copy of the data is necessary
   {
    // copy data, force copying of data CON_Block_t
    bool     temp1 = ownmem;
    unsigned temp2 = size;
    ownmem = true;
    CON_Block_t<T> aux1;
    try {aux1=(*this);} catch (...) {throw;}
    ownmem = temp1;
    // create new object
    if (ownmem)
     {
      clear();
      try {allocate(aux1.layers,aux1.columns,aux1.rows);} catch (...) {clear();throw;}
     }
    else
     {
      T *aux2 = data;
      clear();
      try {assign(aux1.layers,aux1.columns,aux1.rows,aux2);} catch (...) {clear();throw;}
      size = temp2;
     }
    // copy data
    try
     {
      int i,j,k;
      for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (*this)[i][j][k] = aux1[i][k][j];
     }
    catch (...) {clear();throw;}
   }
 }

// transpose CON_Block_t, exchange rows and layers

template<class T> void CON_Block_t<T> :: transpose_rowlay()
 {
  if (rows==layers) // in the symmetric case, we can transpose without actually copying the data
   {
    int i,j,k;
    T   a;
    for (i=0;i<layers;i++) for (j=0;j<i;j++) for (k=0;k<columns;k++) {a=(*this)[i][j][k];(*this)[i][j][k]=(*this)[j][i][k];(*this)[j][i][k]=a;}
   }
  else // for the unsymmetric case an intermediate copy of the data is necessary
   {
    // copy data, force copying of data CON_Block_t
    bool     temp1 = ownmem;
    unsigned temp2 = size;
    ownmem = true;
    CON_Block_t<T> aux1;
    try {aux1=(*this);} catch (...) {throw;}
    ownmem = temp1;
    // create new object
    if (ownmem)
     {
      clear();
      try {allocate(aux1.rows,aux1.layers,aux1.columns);} catch (...) {clear();throw;}
     }
    else
     {
      T *aux2 = data;
      clear();
      try {assign(aux1.rows,aux1.layers,aux1.columns,aux2);} catch (...) {clear();throw;}
      size = temp2;
     }
    // copy data
    try
     {
      int i,j,k;
      for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (*this)[i][j][k] = aux1[j][i][k];
     }
    catch (...) {clear();throw;}
   }
 }

// transpose CON_Block_t, exchange layers and columns

template<class T> void CON_Block_t<T> :: transpose_collay()
 {
  if (layers==columns) // in the symmetric case, we can transpose without actually copying the data
   {
    int i,j,k;
    T   a;
    for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<i;k++) {a=(*this)[i][j][k];(*this)[i][j][k]=(*this)[k][j][i];(*this)[k][j][i]=a;}
   }
  else // for the unsymmetric case an intermediate copy of the data is necessary
   {
    // copy data, force copying of data CON_Block_t
    bool     temp1 = ownmem;
    unsigned temp2 = size;
    ownmem = true;
    CON_Block_t<T> aux1;
    try {aux1=(*this);} catch (...) {throw;}
    ownmem = temp1;
    // create new object
    if (ownmem)
     {
      clear();
      try {allocate(aux1.columns,aux1.rows,aux1.layers);} catch (...) {clear();throw;}
     }
    else
     {
      T *aux2 = data;
      clear();
      try {assign(aux1.columns,aux1.rows,aux1.layers,aux2);} catch (...) {clear();throw;}
      size = temp2;
     }
    // copy data
    try
     {
      int i,j,k;
      for (i=0;i<layers;i++) for (j=0;j<rows;j++) for (k=0;k<columns;k++) (*this)[i][j][k] = aux1[k][j][i];
     }
    catch (...) {clear();throw;}
   }
 }

// mirror CON_Block_t

template<class T>
inline void CON_Block_t<T> :: mirror_horizontal()
 {
  // if CON_Block_t is empty, do nothing
  if (exist())
   {
    int i,j,k;
    T   aux;
    for (k=0;k<layers;k++) for (i=0;i<rows;i++)
     {
      for (j=0;j<columns-j-1;j++)  {aux = (*this)[k][i][j];(*this)[k][i][j]=(*this)[k][i][columns-j-1];(*this)[k][i][columns-j-1]=aux;}
     }
   }
 }

template<class T>
inline void CON_Block_t<T> :: mirror_vertical()
 {
  // if CON_Block_t is empty, do nothing
  if (exist())
   {
    int i,j,k;
    T   aux;
    for (k=0;k<layers;k++) for (i=0;i<rows-i-1;i++)
     {
      for (j=0;j<columns;j++)  {aux = (*this)[k][i][j];(*this)[k][i][j]=(*this)[k][rows-i-1][j];(*this)[k][rows-i-1][j]=aux;}
     }
   }
 }

template<class T>
inline void CON_Block_t<T> :: mirror_depth()
 {
  // if CON_Block_t is empty, do nothing
  if (exist())
   {
    int i,j,k;
    T   aux;
    for (k=0;k<layers-k-1;k++) for (i=0;i<rows;i++)
     {
      for (j=0;j<columns;j++)  {aux = (*this)[k][i][j];(*this)[k][i][j]=(*this)[layers-1-k][i][j];(*this)[layers-1-k][i][j]=aux;}
     }
   }
 }

// concatination of two blocks (horizontal, join columns of each layer)

template<class T>
inline CON_Block_t<T> CON_Block_t<T> :: operator|(const CON_Block_t<T>& org)
 {
  CON_Block_t<T> target;
  if (exist() && org.exist())
   {
    // test size of blocks
    if (org.rows!=rows || org.layers!=layers || !org.data || !data)
        throw Error("CON_Block_t::operator|","Dimensions of operands do not match.",ut_error_array_dimmatch);
    // create target object
    CON_Block_t<T> target;
    try {target.allocate(layers,rows,columns+org.columns);}
    catch (...) {throw;}
    // make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
    // fill in values
    try
     {
      int i,j,k,l;
      for (l=0;l<layers;l++) for (i=0;i<rows;i++)
       {
        for (j=0,k=0;j<columns;j++,k++)   target[l][i][k] = (*this)[l][i][j];
//        for (j=k=0;j<columns;j++,k++)   target[l][i][k] = (*this)[l][i][j];
        for (j=0;j<org.columns;j++,k++) target[l][i][k] = org[l][i][j];
       }
     }
    catch (...) {throw;}
   }
  return target;
 }

// concatination of two blocks (vertical, join rows of each layer)

template<class T>
inline CON_Block_t<T> CON_Block_t<T> :: operator&(const CON_Block_t<T>& org)
 {
  CON_Block_t<T> target;
  if (exist() && org.exist())
   {
    // test size of blocks
    if (org.columns!=columns || org.layers!=layers || !org.data || !data)
        throw Error("CON_Block_t::operator&","Dimensions of operands do not match.",ut_error_array_dimmatch);
    // create target object
    try {target.allocate(layers,rows+org.rows,columns);}
    catch (...) {throw;}
    // make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
    // fill in values
    try
     {
      int i,j,k,l;
      for (l=0;l<layers;l++) for (i=0;i<columns;i++)
       {
        for (j=0,k=0;j<rows;j++,k++)   target[l][k][i] = (*this)[l][j][i];
//        for (j=k=0;j<rows;j++,k++)   target[l][k][i] = (*this)[l][j][i];
        for (j=0;j<org.rows;j++,k++) target[l][k][i] = org[l][j][i];
       }
     }
    catch (...) {throw;}
   }
  return target;
 }

// concatination of two blocks (depth, join layers)

template<class T>
inline CON_Block_t<T> CON_Block_t<T> :: operator^(const CON_Block_t<T>& org)
 {
  CON_Block_t<T> target;
  if (exist() && org.exist())
   {
    // test size of blocks
    if (org.rows!=rows || org.columns!=columns || !org.data || !data)
        throw Error("CON_Block_t::operator^","Dimensions of operands do not match.",ut_error_array_dimmatch);
    // create target object
    try {target.allocate(layers+org.layers,rows,columns);}
    catch (...) {throw;}
    // make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
    // fill in values
    try
     {
      int i,j,k,l;
      for (l=0;l<columns;l++) for (i=0;i<rows;i++)
       {
        for (j=0,k=0;j<layers;j++,k++)   target[k][i][l] = (*this)[j][i][l];
//        for (j=k=0;j<layers;j++,k++)   target[k][i][l] = (*this)[j][i][l];
        for (j=0;j<org.layers;j++,k++) target[k][i][l] = org[j][i][l];
       }
     }
    catch (...) {throw;}
   }
  return target;
 }

// casting operators
template<class T>
inline CON_Block_t<T> :: operator unsigned short() const
{
    CON_Block_t<unsigned short> target;

    try {target.allocate(layers,rows,columns);}
    catch (...) {throw;}

    try
    {
    int i,j,k;
    for(i=0;i<layers;i++)
        for(j=0;j<rows;j++)
            for(k=0;k<columns;k++)
                target[i][j][k] = (*this)[i][j][k];
    }
    catch (...) {throw;}

    return target;
}

// harakiri function: do not use unless authorised
template<class T>
inline T* CON_Block_t<T> :: getdata()
{
    return data;
}

template<class T>
inline const T* CON_Block_t<T> :: getdata() const
{
    return data;
}

//------------------------------------------------------------------------------
// private functions

// initialise unallocated object

template<class T>
inline void CON_Block_t<T> :: init()
 {
  data = NULL;
  datalayers = NULL;
  layers = rows = columns = size = 0;
  ownmem = false;
  reset_ownmem = false;
 }

// common code used in assignment operator and copy constructor

template<class T>
void CON_Block_t<T> :: copy(const CON_Block_t<T>& org)
 {
  clear();
  // copy dimensions and memory flags
  layers = org.layers;
  rows = org.rows;
  columns = org.columns;
  size = org.size;
  if (org.reset_ownmem) ownmem = true;
  else                  ownmem = org.ownmem;
  if (org.ownmem || !org.reset_ownmem) // TRK11022002, if not ownmem=false AND reset_ownmem=true,
                                       // as used in temporary objects that are put on stack by operator functions
   {
    // try to allocate memory, if unsuccessful, throw exception
    if (org.layers*org.rows*org.columns)
     try {data = new T[org.layers*org.rows*org.columns];}
     catch (std::bad_alloc) // seems to be Borland specific (otherwise std::bad_alloc)
     {
      throw Error("CON_Block_t::copy()","Memory could not be allocated.",ut_error_array_alloc);
     }
    if (org.layers)
     try {datalayers = new CON_Matrix_t<T>[org.layers];}
     catch (std::bad_alloc) // seems to be Borland specific (otherwise std::bad_alloc)
     {
      clear();
      throw Error("CON_Block_t::copy()","Memory could not be allocated.",ut_error_array_alloc);
     }
    // copy data
    try
     {
      int ii;
      for (ii=0;ii<org.layers;ii++) datalayers[ii].assign(org.rows,org.columns,data+org.hashtable[ii]);
      memcpy(data,org.data,size*sizeof(T));  //TRK11022002
      ownmem = true;                           //TRK11022002
      // for (ii=0;ii<org.layers*org.rows*org.columns;ii++) data[ii] = org.data[ii];
     }
    catch (...)
     {
      clear();
      throw;
     }
   }
  else // if the original does not own its data
   {
    data = org.data;
    // try to allocate memory, if unsuccessful, throw exception
    if (org.layers)
     try {datalayers = new CON_Matrix_t<T>[org.layers];}
     catch (std::bad_alloc) // seems to be Borland specific (otherwise std::bad_alloc)
     {
      throw Error("CON_Block_t::copy()","Memory could not be allocated.",ut_error_array_alloc);
     }
    // copy data
    try
     {
      for (int ii=0;ii<org.layers;ii++) datalayers[ii].assign(org.rows,org.columns,data+org.hashtable[ii]);
     }
    catch (...)
     {
      clear();
      throw;
     }
   }
  try {hashtable=org.hashtable;}
  catch (...)
   {
    clear();
    throw;
   }
 }

template <class T>
T CON_Block_t<T>::getMaxi() const
{
    T tMaxi = 0;
    if (exist())
    {
        tMaxi = datalayers[0].getMaxi();
        for (long nSample = 1; nSample < layers; ++nSample)
            tMaxi = max(tMaxi, datalayers[nSample].getMaxi());
    }

    return tMaxi;
}

template <class T>
T CON_Block_t<T>::getMini() const
{
    T tMini = 0;
    if (exist())
    {
        tMini = datalayers[0].getMini();
        for (long nSample = 1; nSample < layers; ++nSample)
            tMini = min(tMini, datalayers[nSample].getMini());
    }

    return tMini;
}

/////////////////////////////////////////////////////////////////////////////
// CON_Block_t Serialization
//
#ifdef WIN32
template<class T>
inline void CON_Block_t<T>::Serialize(CON_AbstractPersistenceArchive_i& ar)
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
            throw runtime_error(_T("badClass"));
#endif
        }

        ar << layers;
        ar << rows;
        ar << columns;
        ar << size;
        ASSERT(size == layers * rows * columns);
        ar.Write(data, layers * columns * rows  * sizeof(T));
    }
    else
    {
        ar >> layers;
        ar >> rows;
        ar >> columns;
        ar >> size;
        ASSERT(size == layers * rows * columns);

        allocate(layers, rows, columns);
        reset_ownmem = false;

        ar.Read(data, layers * rows * columns * sizeof(T));
    }
}

template<class T>
inline CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const CON_Block_t<T>& block)
{
    const_cast< CON_Block_t<T>& >(block).Serialize(ar);

    return ar;
}

template<class T>
inline CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, CON_Block_t<T>& block)
{
    block.Serialize(ar);

    return ar;
}
#endif  // WIN32
#endif    // __utBlock_t_H__
