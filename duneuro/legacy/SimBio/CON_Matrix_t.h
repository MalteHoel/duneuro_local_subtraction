//$42   6 july 2007 Steinstraeter  Comments and and i/o methods added. Some bugs fixed.
//$41   15.06.2006  Moritz D.   multiply_matrices_expect_symetric_asresult, virtual void allocate(long,long); for loreta which needs big matrices
//$40   18.02.2005  Markus G.   added dyadicProduct, matrix from the product of two vectors
//$39   15.02.2005  Markus G.   bug in allocate: if old size == new size, dimensions were not set
//$38   18.01.2005  Markus G.   added cholesky inversion (for positive definite matrices)
//$37   10.08.2003  TRK         added mean and sum functions
//$36   31.07.2003    TRK            added exceptions and elementwise operators
//$35   03.01.2003  TRK         new remove functions
//$34    21.03.2002  Anwander A. size was not set in the function assign
//$33    12.03.2002    Anwander A. changed for (j=k=0;j<layers;j++,k++) to for (j=0,k=0;j<layers;j++,k++)
//$32    12.02.2002    TRK            blockcopy with memcpy
//$31    12.02.2002    TRK         copy function makes hardcopy always, if reset_ownmem=false
//$30    12.02.2002    TRK            reset_ownmem = false in assign function
//$29    24.01.2001    A. Anwander    #include <stdexcept> for runtime_error on SGI
//$29    19.02.2001    Frank N.    Added multiply(const CON_Matrix_t<T>& in1, const CON_Vector_t<T>& in2, CON_Vector_t<T>& out)
//$28    13.11.2000    Frank N.    Made some functions virtual
//$27    07.11.2000    Frank Z.    changed copy() from private to protected, made allocate() virtual
//$26    14.08.2000    Frank N.    Changed return value for const CON_Vector_t<T>& operator[](int) const
//$25   28.07.2000  Matthias D. Made changes on multiplications with a diagonalmatrix, diagonal matrix values in vector
//$24   25.07.2000  Matthias D. Added multiplications with a diagonalmatrix, diagonal matrix values in vector
//$23    10.07.2000    Matthias D.    Removed Capital letters from file names
//$22    05.07.2000    Frank N.    Adapted for Simbio
//$21    11.04.2000    Frank Z.    Added getMini and getMaxi()
//$20    07.04.2000    Frank N.    Added multiply()
//$19    06.04.2000    Frank N.    Added removerow() and removecolumn()
//$18    26.03.2000    Frank N.    Smart matrix mult.
//$17    23.03.2000    Frank N.    Two operator*()'s became const
//$16    15.03.2000    Frank N.    Added gauss(), swaprows() and getMaxPosition()
//$15    14.03.2000    Frank N.    Handle the NOT "ownmem" case in the Serialize() function
//$14    13.03.2000    Frank N.    Worked on performance improvements
//$13    24.02.2000    Frank N.    Changed define for bad_allcoc
//$12    21.02.2000    Frank Z.    bug fix in allocation routine
//$11    17.02.2000    Frank N.    Derived from CON_PersistentClass_t
//$10    20.09.1999    Frank Z.    bug fix in concatenate operator &, |
//$9    27.08.1999    Frank Z.    Added #include <CON_Vector_t.h>
//$8    26.08.1999    Frank Z.    implemented column
//$7    20.08.1999    Frank N.    Fixed setcolumn()
//$6    18.08.1999    Frank Z.    implemented setrow, setcolumn
//$5    09.08.1999    TRK            Changed to #include <stdexcept>
//$4    01.07.1999    Frank N.    Serialization is now based on CON_AbstractPersistenceArchive_i
//$3    21.06.1999    Frank Z.    bug fixed
//$2    10.05.1999    Frank N.    added serialization - for now for CArchive only
//$1    05.05.1999    Frank N.    split from ut_array_t.h
//-------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////
//
//    NeuroFEM license:
//    =================
//    Copyright (c) 2007 by
//    Dr.Carsten Wolters, Dr.Alfred Anwander, Dr.Matthias Duempelmann,
//    Dr.Thomas Knoesche, Moritz Dannhaur, Dr. Olaf Steinstraeter.
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
// File name : CON_Matrix_t.h
//-------------------------------------------------------------------------

#ifndef __CON_MATRIX_T_H__
#define __CON_MATRIX_T_H__

#include <stdexcept>

#include "CON_UtilitiesDef.h"

#include "Error.h"
#include "CON_Vector_t.h"

#include "CON_AbstractPersistenceArchive_i.h"
#include "CON_PersistentClass_t.h"

// added by Steinstraeter: BGN
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// added by Steinstraeter: END

#ifdef WIN32
// Forward declarations
template <class T>
class std::valarray;
#endif

// Two-dimensional array (CON_Matrix_t)
// This object can be used like an ordinary array, it is however saveguarded
// against out-of-range errors. Additional functions of e.g. CON_Matrix_t algebra are
// available
// comment added by Steinstraeter:
//         (1) For an CON_Matrix_t<T>, the array of T values are initialised
//             by memset(...,0, ...). Therefore, it must be possible to
//             initialize instances of T by sizeof(T) bytes with value 0.
//             Class instances are normally destroyed by this mechanism.
//         (2) The internal block of data (getdata()) may be larger than the
//             actually represented matrix. The method assign_reduced() assigns
//             a block of data but removes some entries from the actually
//             represented matrix. Therefore getdata(), interpreted as a continuous
//             block of data is different from the represented matrix. Especially
//             height() and length() does not describe the data returned by getdata().
//             To use getdata() as an interface to FORTRAN routines (LAPACK) is therefore
//             dangerous.
//------------------------------------------------------------------------------
// The following functions are available. X,Y,Z are matrices, i,j,k,rows,cols
// integers, a,b,c double/floats.
//
// description                        example                    comment
// ------------                       ---------------------      ---------------
// creation                           CON_Matrix_t<float> X;
// creation as copy                   CON_Matrix_t<long> X(Y);
// creation and allocation            CON_Matrix_t<char> X(rows,cols);
// allocation                         X.allocate(rows,cols)
//    comment added by Steinstraeter: X.allocate(0,0) is not the same as X.clear(), instead, X is unchanged after X.allocate(0, 0).
// assignment of data                 int i[22];X.assign(2,11,i);
// deallocation                       clear(X);
// test of existence                  if (exist(X))             // true if CON_Matrix_t is allocated with both dimensions > 0
// dimensions of CON_Matrix_t               i = length(X)*height(X);  // length is number of columns, height number of rows
// reference to rows and elements     CON_Vector_t<int> z = X[3];
//                                    int a = X[3][1];
// reference to columns               CON_Vector_t<int> x = X.column(2)
// overwrite row CON_Vector_t               X.setrow(x,3);             // the length of the CON_Vector_t and the number of columns / rows
// overwrite column CON_Vector_t            X.setcolumn(x,4);          // of the CON_Matrix_t must match, the given indices must be within range
// assignment                         CON_Matrix_t<char> X = Y;
// addition (element-wise)            X = Y + Z; Z += Y;         // operands must have equal dimensions
// subtraction (element-wise)         X = Y - Z; Z -= Y;         // operands must have equal dimensions
// matrix multiplication              CON_Matrix_t<int> X(2,3),Y(3,5); // columns of first and rows of second operand must be equal
// element-wise multiplication          X = A%B
// element-wise division              X = A/B
//                                    CON_Matrix_t<int> Z=X*Y;         // new CON_Matrix_t has dimensions (2,5)
// multiplication with scalar         X *= a; Y = Z * b;         // NOT Y = b * Z !!!
// equality / unequality              if (X == Y && Y != Z)
// Frobenius norm (yields double)     double a = norm(X);
//    comment added by Steinstraeter: WARNING: norm() returns the square of the Frobenius norm not the Frobenius norm !!!!!!
// mean over all elements              double a = mean(X);
// sum of all elements                  double a = sum(X);
// transpose                          transpose(X);              // NO assignment (instead of Y=transpose(X), use Y=X;transpose(Y)),
//                                                               // CAUTION: CON_Matrix_t is manipulated
//                                                               // CAUTION: in asymmetric case the CON_Matrix_t is temporarily copied
// mirror                             mirror_horizontal(X);      // reversing column indices
//                                    mirror_vertical(Z);        // reversing row indices
//                                                               // CAUTION: as with transpose, no assignment is possible. This
//                                                               //          is a safety measure, because the CON_Matrix_t is manipulated
// concatination                      X = (Y | Z) & W            // first concatinates columns of Y and Z, and then concatinates the rows of the result
//                                                               // with the rows of W, dimensions must match.
// dyadic product                     dyadicProduct(v1, v2);     // A = v1 * v2^T, create matrix from multiplication of two vectors
// Cholesky-inverse                   invertCholesky();          // calculate inverse (^-1) of a symmetric positive definite matrix
//



// added by Steinstraeter: BGN
// ========================================================================================================================
// Some Extensions
// ===============
//
// Legend:
// -------
//
//    CON_Matrix_t<T> M, S;
//    ostream out;
//
// Class information:
// ------------------
//
//    char info = M.getElementType();
//      Information about the template parameter T used for M.
//      Class method.
//      info = 'd' if T = double, info = 'f' if T = float, info = '?' otherwise.
//
// Data Access:
// ------------
//
//    CON_Vector_t<T> vec = M(i);
//      =  M[i], but without range check as long as RANGECHECKDEBUG is not set
//    T value = M(i)(j);
//      = M[i][j], but without range check as long as RANGECHECKDEBUG is not set
//    T value = M(i, j);
//      = M(i)(j)
//
// Resize:
// -------
//    M.allocate_or_clear(new_rows, new_columns)
//      calls clear() if new_columns <= 0 or new_columns <= 0
//      and allocate(new_rows, new_columns) otherwise.
//
// Data Exchange:
// --------------
//
//    M.swap(S)
//    Exchange content of M and S.
//
// I/O methods:
// ------------
//
//    out << M;
//      writes matrix M to out (class ostream, e.g. cout or cerr), matlab notation
//
//    std::string s = M.toString(maxLen, setw_value, format, setprecision_value, use_cr);
//    std::string s = M.toString(maxLen, setw_value, format, setprecision_value);
//    std::string s = M.toString(maxLen, setw_value, format);
//    std::string s = M.toString(maxLen, setw_value);
//    std::string s = M.toString(maxLen);
//    std::string s = M.toString();
//      converts M into string, matlab notation.
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
//         use_cr:               rows are separated by '\n' if use_cr == true,
//                               default: false.
//
//    bool success = M.save(filename);
//      writes M to binary file filename.
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
//        <rows> = M.height(), <columns> = M.length(),
//        <type> = 'f' if T = float, 'd' otherwise.
//      A corresponding matlab function utMatrixLoad.m exists in <SimBioHome>/utilities/matlab.
//      success = false indicates error.
//    bool success = M.load(filename);
//      read M from binary file filename.
//      Format: see M.save(filename).
//              <rows> -> M.height(), <columns> -> M.length(),
//      A corresponding matlab function utMatrixSave.m exists in <SimBioHome>/utilities/matlab.
//      success = false indicates error (in this case: M is empty (M.exist() = false)).
//
// BLAS / LAPACK convenience functions:
// ------------------------------------
//
//    rebuildLowerTriang(M);
//    For symmetric matrices some BLAS / LAPACK function only access the
//    upper or lower part of the passed matrix. If M is a result of such
//    an operation with UPLO = 'L' (lower triangle), only the values of the
//    upper triangle (implicit transpose after interpreting as CON_Matrix_t)
//    are valid. rebuildLowerTriang(M) copys the upper triangle to the lower
//    one.
//
//    packSymMatrix(unpackedMatrix, packedVector)
//      Compact storage of symmetric matrices in a vector: pack operation
//      The upper triangle of unpackedMatrix is packed by rows:
//         unpackedMatrix:   NxN
//                           a11 a12 a13 a14
//                            *  a22 a23 a24
//                            *   *  a33 a34
//                            *   *   *  a44
//         packedVector:     1 x N*(N+1)/2
//                           a11 a12 a13 a14  a22 a23 a24  a33 a34  a44
//      BLAS and LAPACK use the following packing sheme if UPLO = 'L'
//        (lower triangle) is specified, packing by rows:
//       a11  *   *   *
//       a21 a22  *   *   -> a11 a21 a31 a41  a22 a32 a42  a33 a43  a44
//       a31 a32 a33  *
//       a41 a42 a43 a44
//      This leads to an identical result for the packed matrix, if the matrix is symmetric.
//    unpackSymMatrix(packedVector, unpackedMatrix)
//      Compact storage of symmetric matrices in a vector: unpack operation
//      Inverse of packSymMatrix().
//      unpackSymMatrix() only accesses the upper triangle of the matrix unpackedMatrix
//      all other elements are untouched (0 if unpackedMatrix.allocate() changes the size
//      of unpackedMatrix). The remaining elements can be set by rebuildLowerTriang().
//      With n = packedVector.length(), unpackedMatrix becomes a NxN matrix with
//      N = sqrt(2*n + 1/4)-1/2. To avoid this calculation, N can explicitly be specified.
//      N < 0  implicitly uses sqrt(2*n + 1/4)-1/2.
// ========================================================================================================================
// added by Steinstraeter: END



template <class T>
class CON_Matrix_t : public CON_PersistentClass_t< CON_Matrix_t< T > >
{
public:
    CON_Matrix_t();                           // create CON_Matrix_t with 0x0 elements
    CON_Matrix_t(long,long,T=0);                 // create empty CON_Matrix_t of given dimensions
    virtual void allocate(long,long);             // through away and reallocate CON_Matrix_t
  void allocate_or_clear(long new_rows, long new_columns); // added by Steinstraeter

    CON_Matrix_t(const CON_Matrix_t<T>&);           // copy CON_Matrix_t
    virtual ~CON_Matrix_t();                          // delete CON_Matrix_t

    void clear();                       // clear CON_Matrix_t
  void swap(CON_Matrix_t<T> &other);    // added by Steinstraeter: exchange the content of *this and other
    void assign(int,int,T*);            // assign data to CON_Matrix_t
    virtual void assign_reduced(const CON_Matrix_t<T>&,const CON_Vector_t<bool>&,unsigned,unsigned);
                                       // assign other CON_Matrix_t to this CON_Matrix_t with passed selection
    bool exist() const;                 // return allocation status
    unsigned length() const;            // horizontal dimension
    unsigned height() const;            // vertical dimension
    unsigned size_of_data() const;      // added by Steinstraeter: size of the data block accessible by getdata(),
                                      //                         WARNING: If size_of_data() != height()*length() do not use getdata().
  int getMaxPosition(int col) const; // comment added by Steinstraeter: unusable, compiler error.
    CON_Vector_t<T>& operator[](unsigned int);         // reference a CON_Matrix_t row, can be further de-referenced using the CON_Vector_t operator
    const CON_Vector_t<T>& operator[](unsigned int) const;
    CON_Vector_t<T>& operator()(int); // added by Steinstraeter: identical to operator[] but without range check (as long as RANGECHECKDEBUG is not set).
    const CON_Vector_t<T>& operator()(int) const; // added by Steinstraeter: identical to operator[] but without range check (as long as RANGECHECKDEBUG is not set).
  T &operator()(int, int); // added by Steinstraeter: identical to [][] but without range check (as long as RANGECHECKDEBUG is not set).
  const T &operator()(int, int) const; // added by Steinstraeter: identical to [][] but without range check (as long as RANGECHECKDEBUG is not set).
    CON_Vector_t<T> column(int) const;       // get a copy of column
    void setrow(const CON_Vector_t<T>&,int);       // copy CON_Vector_t into ith row of CON_Matrix_t
    void setcolumn(const CON_Vector_t<T>&,int);    // copy CON_Vector_t into ith column of CON_Matrix_t
    void removerow(int);                        // remove the ith row of CON_Matrix_t
    void removecolumn(int);                        // remove the ith column of CON_Matrix_t
    void removerows(int,int);                    // remove from ith row of CON_Matrix_t n rows
    void removecolumns(int,int);                // remove from ith column of CON_Matrix_t n columns
    CON_Matrix_t<T>& operator=(const CON_Matrix_t<T>&);   // assign CON_Matrix_t to another CON_Matrix_t
    CON_Matrix_t<T>& operator=(const T&);           // fill with one value
    CON_Matrix_t<T>& operator+=(const CON_Matrix_t<T>&);  // add another CON_Matrix_t (same size)
    CON_Matrix_t<T>  operator+(const CON_Matrix_t<T>&);
    CON_Matrix_t<T>& operator-=(const CON_Matrix_t<T>&);  // subtract another CON_Matrix_t (same size)
    CON_Matrix_t<T>  operator-(const CON_Matrix_t<T>&);
    CON_Matrix_t<T>  operator*(const CON_Matrix_t<T>&) const;   // multiplication of two matrices (special dimension requirements)
    CON_Matrix_t<T> operator%(const CON_Matrix_t<T>&) const;
    CON_Matrix_t<T> operator/(const CON_Matrix_t<T>&) const;
    static void multiply(const CON_Matrix_t<T>& in1, const CON_Matrix_t<T>& in2, CON_Matrix_t<T>& out); // multiplication of two matrices (special dimension requirements)
    void  multiplyrightdiagonal(const CON_Vector_t<T>& in,CON_Matrix_t<T>& out); // multiplications with a diagonalmatrix on the right side, diagonal matrix values in vector in2
    void  multiplyleftdiagonal(const CON_Vector_t<T>& in, CON_Matrix_t<T>& out); // multiplications with a diagonalmatrix on the left side, diagonal matrix values in vector in1
    CON_Vector_t<T>  operator*(const CON_Vector_t<T>&) const;   // multiplication of the CON_Matrix_t with CON_Vector_t
    static void multiply(const CON_Matrix_t<T>& in1, const CON_Vector_t<T>& in2, CON_Vector_t<T>& out); // multiplication of matrix with vector (special dimension requirements)
    CON_Matrix_t<T>& operator*=(const T&);    // multiply with scalar
    CON_Matrix_t<T>  operator*(const T&)    const;
    bool operator==(const CON_Matrix_t<T>&) const;  // equality
    bool operator!=(const CON_Matrix_t<T>&) const;  // unequality
    double norm() const;                // frobenius norm
    double sum() const;
    double mean() const;
        void multiply_matrices_expect_symetricmatrix_asresult(CON_Matrix_t<T>&in1 ,CON_Matrix_t<T>& out);
     bool invertCholesky();            // invert matrix using Cholesky decomposition
        void setValue(long x, long y,double value);
        double getValue(long x, long y);
    void dyadicProduct(const CON_Vector_t<T> &v1, const CON_Vector_t<T> &v2);
// void invert();                       // compute inverse
#ifdef WIN32
    bool gauss(std::valarray<long double> & x, const std::valarray<long double> & b) const;
#endif
    void transpose();                   // swap rows and columns
    void mirror_horizontal();           // X[i][j] becomes X[i][n-j]
    void mirror_vertical();             // X[i][j] becomes X[m-i][j]
    void swaprows(int row1, int row2);
    CON_Matrix_t<T>  operator|(const CON_Matrix_t<T>&) const;   // concatenate horizontally
    CON_Matrix_t<T>  operator&(const CON_Matrix_t<T>&) const;   // concatenate vertically
    T* getdata();                                            // return pointer to data
                                    // comment added by Steinstraeter: WARNING: If size_of_data() != height()*length() do not use getdata().
    const T* getdata() const;                    // return pointer to data
                                    // comment added by Steinstraeter: WARNING: If size_of_data() != height()*length() do not use getdata().

    T getMaxi() const;
    T getMini() const;

  static char getElementType(); // added by Steinstraeter: 'f' if T = float, 'd' if T = double, '?' otherwise.
  std::string toString(int maxLen = -1,
                       int setw_value = 0,
                       char format = 'g',
                       int setprecision_value = -1,
                       bool use_cr = false) const; // added by Steinstraeter: string representation, matlab notation.
  bool save(const char *filename) const; // added by Steinstraeter: write matrix elements to binary file.
  bool load(const char *filename); // added by Steinstraeter: read matrix elements from binary file.

    // Serialization
#ifdef WIN32
    virtual void Serialize(class CON_AbstractPersistenceArchive_i& ar);
    friend CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const CON_Matrix_t<T>& matrix);
    friend CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, CON_Matrix_t<T>& matrix);
#endif

protected:
   virtual void copy(const CON_Matrix_t<T>&);        // common code used in assignment operator and copy constructor
   void init();                                    // initialise unallocated object

protected:
   unsigned                rows,columns;    // dimensions of CON_Matrix_t
   unsigned                size;           // size of data allocated CON_Block_t (0 if no allocation)
                                  // comment added by Steinstraeter: Not equal to rows * columns. The size of
                                  //                                 data may be larger than the respresented matrix.
                                  //                                 rows, columns, datarows, and hashtable define the
                                  //                                 matrix represented by the CON_Matrix_t object. This
                                  //                                 matrix may be a subset of data.
   CON_Vector_t<T>*        datarows;       // CON_Vector_t objects for rows
   T*                    data;           // data CON_Matrix_t
                              // comment added by Steinstraeter: Data may be larger than the matrix represented by
                              //                                 the CON_Matrix_t object, see above (size).
   bool                    ownmem;         // indicates, whether the allocated memory is owned by this object and
                                        // should be copied and deleted by it
   bool                    reset_ownmem;   // if set, the copy constructor resets ownmem to true
                                // comment added by Steinstraeter:
                                //   If src.reset_ownmem = true and src.ownmem = false, the next dest.copy(src)
                                //   will transfer the responsibility of freeing data to dest.
                                //   This is useful for a temporary vector that should be copied to a function caller
                                //   by return (e.g. operator+()).
   CON_Vector_t<unsigned> hashtable;      // indicates offsets of the row pointers

   DECLARE_PERSISTENT_CLASS( CON_Matrix_t< T > )
};

// some non-element equivalents (just for more clear notation, e.g. if (norm(x)<a) instead of if (x.norm()<a)

template<class T>
inline void clear (CON_Matrix_t<T>& arg)
{
    arg.clear();
}

template<class T>
inline bool exist (CON_Matrix_t<T>& arg)
{
    return arg.exist();
}

template<class T>
inline unsigned length(CON_Matrix_t<T>& arg)
{
    return arg.length();
}

template<class T>
inline unsigned height(CON_Matrix_t<T>& arg)
{
    return arg.height();
}

template<class T>
inline double norm(CON_Matrix_t<T>& arg)
{
    return arg.norm();
}

template<class T>
inline double sum(CON_Matrix_t<T>& arg)
{
    return arg.sum();
}

template<class T>
inline double mean(CON_Matrix_t<T>& arg)
{
    return arg.mean();
}

// added by Steinstraeter:
//   formated stream output, matlab notation
template<class T>
std::ostream &operator<<(std::ostream &s, const CON_Matrix_t<T> &M)
{
  s << "[";
  if (M.exist())
    for (unsigned int i = 0; i < M.height(); i++)
      {
        for (unsigned int j = 0; j < M.length()-1; j++)
          s << M[i][j] << " ";
        s << M[i][M.length()-1];
        if (i < M.height()-1)
          s << "; ";
      }
  s << "]";

  return s;
}

// added by Steinstraeter:
//   Rebuild the lower triangular part (without the diagonal) of the symmetric matrix M.
template<class T>
inline void rebuildLowerTriang(CON_Matrix_t<T> &M)
{
  for (int i = 1; i < M.height(); i++)
    for (int j = 0; j < i; j++)
      M(i,j) = M(j,i);
}

// added by Steinstraeter: Compact storage of symmetric matrices in a vector: pack operation
template<class T>
inline void packSymMatrix(const CON_Matrix_t<T> &unpacked, CON_Vector_t<T> &packed)
{
  int N = min(unpacked.height(), unpacked.length());
  int packedSize = (N*(N+1))/2;
  if (packed.length() != packedSize)
    packed.allocate_or_clear(packedSize);

  int bytes_to_copy = N * sizeof(T);
  char *dest = reinterpret_cast<char*>(packed.getdata());
  for (int i = 0; i < N; i++)
    {
       memcpy(dest, unpacked(i).getdata() + i, bytes_to_copy);
       dest += bytes_to_copy;
       bytes_to_copy -= sizeof(T);
    }
}

// added by Steinstraeter: Compact storage of symmetric matrices in a vector: unpack operation
template<class T>
inline void unpackSymMatrix(const CON_Vector_t<T> &packed, CON_Matrix_t<T> &unpacked, int N = -1)
{
  if (N < 0)
    N = static_cast<int>(sqrt(2*static_cast<double>(packed.length()) + 0.25)); // round(sqrt(2*packed.length() + 0.25) - 0.5)

  if (unpacked.height() != N || unpacked.length() != N)
    unpacked.allocate_or_clear(N, N);

  int bytes_to_copy = N * sizeof(T);
  char *src = reinterpret_cast<char*>(const_cast<CON_Vector_t<T>&>(packed).getdata());
  for (int i = 0; i < N; i++)
    {
       memcpy(unpacked(i).getdata() + i, src, bytes_to_copy);
       src += bytes_to_copy;
       bytes_to_copy -= sizeof(T);
    }
}

// template<class T> inline void     invert (CON_Matrix_t<T>& arg) {arg.invert();}

template<class T>
inline void transpose (CON_Matrix_t<T>& arg)
{
    arg.transpose();
}

template<class T>
inline void mirror_horizontal(CON_Matrix_t<T>& arg)
{
    arg.mirror_horizontal();
}

template<class T>
inline void mirror_vertical (CON_Matrix_t<T>& arg)
{
    arg.mirror_vertical();
}

// constructor: just initialise

template<class T>
CON_Matrix_t<T> :: CON_Matrix_t()
{
    init();
}

// constructor and allocator: allocate CON_Matrix_t with i x j elements of type T

template<class T>
CON_Matrix_t<T> :: CON_Matrix_t(long i,long j,T a)
{
    init();
    allocate(i,j);
    (*this) = a;
}

template<class T>
inline void CON_Matrix_t<T> :: allocate(long i,long j)
{
// Check for positive value first
    ASSERT(i >= 0 || j >= 0);
    if (i<0 || j<0)
        throw Error("CON_Matrix_t::allocate()","Wrong index!",ut_error_array_index);
    if (i == 0 || j == 0)
        return;

// Check for same dimensions afterwards, if so just reuse the old memory
    if (i == height() && j == length() && ownmem == true && data != NULL)
    {
        ::ZeroMemory(data, sizeof(T) * size);
        return;
    }

// Allocate data CON_Block_t
    clear();

    if (i > 0 && j > 0)
    {
        try
        {
            data = new T[i*j];
            datarows = new CON_Vector_t<T>[i];
            hashtable.allocate(i);

            for (long ii = 0; ii < i; ii++)
            {
                datarows[ii].assign(j, data + ii * j);
                hashtable[ii] = ii * j;
            }
            ownmem = true;
            ::ZeroMemory(data, sizeof(T) * i * j);
        }
        // in case of failure
        catch (std::bad_alloc&)
        {
            clear();
            ASSERT(FALSE);
            throw Error("utMatrix::allocate()","Allocation failed!",ut_error_array_alloc);
        }
    }

    rows = i;
    columns = j;
    size = i * j;
}

// added by Steinstraeter: clear() if new_rows or new_columns <= 0, otherwise allocate().
template<class T>
inline void CON_Matrix_t<T>::allocate_or_clear(long new_rows, long new_columns)
{
 if (new_rows <= 0 || new_columns <= 0)
   clear();
 else
   allocate(new_rows, new_columns);
}


// contructor: copy another CON_Matrix_t
template<class T>
inline double CON_Matrix_t<T> :: getValue(long x, long y)
{
 return (*this)[x][y];
}

template<class T>
inline void CON_Matrix_t<T> :: setValue(long x, long y, double value)
{
 (*this)[x][y]=value;
}

template<class T>
inline void CON_Matrix_t<T> :: dyadicProduct(const CON_Vector_t<T> &v1,
    const CON_Vector_t<T> &v2)
{
    int l1=v1.length();
    int l2=v2.length();
    allocate(l1, l2);
    for (int i=0; i<l1; i++)
        for (int j=0; j<l2; j++)
            datarows[i][j] = v1[i] * v2[j];
}

// simple matrix multiplication with expected symetric structure as result
// comment added by Steinstraeter: *this and in1 must be square matrices with identical size.
template<class T>
inline void CON_Matrix_t<T>:: multiply_matrices_expect_symetricmatrix_asresult(CON_Matrix_t<T>&in1 ,CON_Matrix_t<T>& out)
{

 long len=(* this ).length();
 if(len != in1.height())
 {
   return ;
 }
 CON_Vector_t< double> vec1,vec2;

 double tmp;

 out.allocate(len,len);


 for(long i=0;i<len;i++)
 {
  vec1=in1[i];
  for(long j=0;j<len;j++)
   {
   if(j>=i)
    {
     vec2=(* this).column(j);
     tmp=vec1*vec2;
     out[i][j]=tmp;
     if (j!=i) out[j][i]=tmp;
    }
   }

 }

}

 //invert matrix, using Cholesky Decomposition
//only works for symmetric positive definite matrices
//no check for symmetry, the result is meaningless if matrix was not symmetric
//returns false, if matrix isn't square or not positive definite, else true
//the algorithms used are based on [NR, Numerical Recipes in C++, 2nd edition],
//chapter 2.9, pages 99-101
template<class T>
bool CON_Matrix_t<T> :: invertCholesky()
{
    if (!exist() || height() != length())        //check if matrix is square
        return false;

    //a symmetry-check is not performed, to save time

    int n = length();                            //n == height == length

    //Cholesky decomposition of A
    //calculate lower triangular matrix L with A = L * L^T (L^T is transpose of L)
    //Algorithm from [NR::choldc] p. 100
    for (int i=0; i<n; i++)
    {
        for (int j=i; j<n; j++)
        {
            T sum = data[i*n + j];
            for (int k = i-1; k>=0; k--)
                sum -= data[i*n + k] * data[j*n + k];

            if (i==j)
            {
                if (sum <= 0)
                    return false;                //matrix is not positive definite
                data[i*n + i] = sqrt(sum);
            }
            else
                data[j*n + i] = sum / data[i*n + i];
        }
    }


     //calculate L^-1 (the inverse of L)
    //Algorithm from [NR] p.101, but elements are also stored in the upper
    //triangle for better performance here and in the following multiplication
    //Like this you can multiply row-wise, resulting in faster memory-access

    for (int i=0; i<n; i++)
    {
        data[i*n + i] = 1.0 / data[i*n + i];
        for (int j = i+1; j<n; j++)
        {
            T sum = 0;
            for (int k=i; k<j; k++)
                sum -= data[j*n + k] * data[i*n + k];    //multiply row-wise

            T result = sum / data[j*n + j];
            data[j*n + i] = result;                //store in lower triangle
            data[i*n + j] = result;                //store in upper triangle
        }
    }

    //calculate A^-1 (the inverse of the original matrix) by A^-1 = L-1^T * L^-1
    for (int i=0; i<n; i++)
    {
        for (int j=i; j<n; j++)
        {
            T sum = 0;
            for (int k=j; k<n; k++)
                sum += data[i*n + k] * data[j*n + k];
            data[i*n + j] = sum;
            data[j*n + i] = sum;
        }
    }

    return true;
}


template<class T>
CON_Matrix_t<T> :: CON_Matrix_t(const CON_Matrix_t<T>& org)
{
    init();
    copy(org);
}

// destructor and clear function: destroy CON_Matrix_t and free ressources

template<class T>
CON_Matrix_t<T> :: ~CON_Matrix_t()
{
    clear();
}

template<class T>
inline void CON_Matrix_t<T> :: clear()
{
// in case of failure, dont do anything
    if (datarows)
    {
        try
        {
            delete [] datarows;
        }
        catch(...)
        {
            ASSERT(FALSE);
            throw Error("CON_Matrix_t::clear()","Freeing of memory failed!",ut_error_array_alloc);
        };
    }
    if (ownmem && data)
    try
    {
        delete [] data;
    }
    catch(...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::clear()","Freeing of memory failed!",ut_error_array_alloc);
    };

    init();
}

// assign array of type T to CON_Matrix_t

template<class T>
inline void CON_Matrix_t<T> :: assign(int i,int j, T* orgdata)
{
    if (i<0 || j<0)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::assign(...)","Wrong index!",ut_error_array_index);
        return;
    }

    if (i==0 || j==0)
        return;

    clear();
    data = orgdata;
    try
    {
        datarows = new CON_Vector_t<T>[i];
        hashtable.allocate(i);
        for (int ii=0;ii<i;ii++)
        {
            datarows[ii].assign(j,data+ii*j);
            hashtable[ii] = ii*j;
        }
    }
    catch(...)
    {
        clear();
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::assign(...)","Assignment failed!",ut_error_unknown);
    };
    rows = i;
    columns = j;
    size = i * j;          // Anwander A. 21/03/2002 size was not set
    ownmem = false;
    reset_ownmem = false; // TRK11022002, should not own the memory later either
}

// assign array of type T to CON_Matrix_t with passed hashtable

template<class T>
void CON_Matrix_t<T> :: assign_reduced(const CON_Matrix_t<T>& full,const CON_Vector_t<bool>& select,unsigned offset,unsigned samples)
{
// if full CON_Matrix_t is empty, just return
    if (!full.exist())
        return;

// if selection CON_Vector_t is empty, just copy full
    if (!select.exist())
    {
        clear();
        copy(full);
        return;
    }

// if number of rows of full CON_Matrix_t and length of selection is different
    if (select.length()!=full.height())
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::assign_reduced(...)","Wrong selection!",ut_error_array_selectrows);
        return;
    }

// if offset or length disagree with dimensions of full CON_Matrix_t
    if (offset < 0 || offset > full.length()-1 || samples < 1 || samples > full.length()-offset)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::assign_reduced(...)","Wrong offset or length!",ut_error_array_selectrows);
        return;
    }

    clear();
// determine numbers of rows and columns
    unsigned int i,j;
    rows=0,columns=samples;
    for (i=0;i<full.height();i++)
        if (select[i])
            rows++;

    if (!rows || !columns)
    {
        clear();
        return;
    }

// create new object
    size = full.size;
    data = full.data;
    try
    {
        datarows = new CON_Vector_t<T>[rows];
        hashtable.allocate(rows);
        for (i=j=0;i<full.height();i++)
        {
            if (select[i])
            {
                datarows[j].assign(columns,data+full.hashtable[i]+offset);
                hashtable[j] = full.hashtable[i] + offset;
                j++;
            }
        }
    }
    catch (std::bad_alloc&)
    {
        clear();
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::assign_reduced(...)","Allocation failed!",ut_error_array_alloc);
    }

    ownmem = false;
    reset_ownmem = false; // TRK11022002, should not own the memory later either
}

// return allocation status

template<class T>
inline bool CON_Matrix_t<T> :: exist() const
{
    if (rows && columns && data && datarows)
        return true;
    else
        return false;
}

// return dimensions

template<class T>
inline unsigned CON_Matrix_t<T> :: length() const
{
    return columns;
}

template<class T>
inline unsigned CON_Matrix_t<T> :: height() const
{
    return rows;
}

template<class T>
inline unsigned CON_Matrix_t<T> :: size_of_data() const
{
    return size;
}

// comment added by Steinstraeter: unusable, compiler error.
template<class T>
inline int CON_Matrix_t<T> :: getMaxPosition(int col) const
{
    long double valmax;
    int maxvalpos = col;

    valmax = Abs((*this)[col][col]);
    for (int i = col+1; i < height(); i++)
    {
        if (Abs((*this)[i][col]) > valmax)
        {
            valmax = Abs((*this)[i][col]);
            maxvalpos = i;
        }
    }

    return maxvalpos;
}

// overload brackets to reference the ith row of the CON_Matrix_t

template<class T>
inline CON_Vector_t<T>& CON_Matrix_t<T> :: operator[](unsigned int i)
{
    static CON_Vector_t<T> errorVar;

    if (i<0 || i>=rows)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator[]","Wrong index!",ut_error_array_index);
    }
    return datarows[i];
}

template<class T>
inline const CON_Vector_t<T>& CON_Matrix_t<T> :: operator[](unsigned int i) const
{
    static CON_Vector_t<T> errorVar;

    if (i<0 || i>=rows)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator[]","Wrong index!",ut_error_array_index);
    }
    return datarows[i];
}

// added by Steinstraeter: exchange the content of *this and other
template<class T>
void CON_Matrix_t<T>::swap(CON_Matrix_t<T> &other)
{
  std::swap(rows, other.rows);
  std::swap(columns, other.columns);
  std::swap(size, other.size);
  std::swap(datarows, other.datarows);
  std::swap(data, other.data);
  std::swap(ownmem, other.ownmem);
  std::swap(reset_ownmem, other.reset_ownmem);
  hashtable.swap(other.hashtable);
}

// added by Steinstraeter: identical to operator[] but without range check as long as RANGECHECKDEBUG is not set
template<class T>
inline CON_Vector_t<T>& CON_Matrix_t<T> :: operator()(int i)
{
#ifdef RANGECHECKDEBUG
    if (i<0 || i>=rows)
    {
          std::cerr << "CON_Matrix_t::operator(): " << "Wrong index!" << std::endl;
      exit(1);
    }
#endif
    return datarows[i];
}

// added by Steinstraeter: identical to operator[] but without range check as long as RANGECHECKDEBUG is not set
template<class T>
inline const CON_Vector_t<T>& CON_Matrix_t<T> :: operator()(int i) const
{
#ifdef RANGECHECKDEBUG
    if (i<0 || i>=rows)
    {
          std::cerr << "CON_Matrix_t::operator(): " << "Wrong index!" << std::endl;
      exit(1);
    }
#endif
    return datarows[i];
}

// added by Steinstraeter: identical to [][] but without range check as long as RANGECHECKDEBUG is not set
template<class T>
inline T &CON_Matrix_t<T>::operator()(int i, int j)
{
#ifdef RANGECHECKDEBUG
    if (i<0 || i>=rows)
    {
          std::cerr << "CON_Matrix_t::operator(): " << "Wrong index!" << std::endl;
      exit(1);
    }
#endif
    return datarows[i](j);
}

// added by Steinstraeter: identical to [][] but without range check as long as RANGECHECKDEBUG is not set
template<class T>
inline const T &CON_Matrix_t<T>::operator()(int i, int j) const
{
#ifdef RANGECHECKDEBUG
    if (i<0 || i>=rows)
    {
          std::cerr << "CON_Matrix_t::operator(): " << "Wrong index!" << std::endl;
      exit(1);
    }
#endif
    return datarows[i](j);
}

// overload assignment operator

template<class T>
inline CON_Matrix_t<T>& CON_Matrix_t<T> :: operator=(const CON_Matrix_t<T>& org)
{
    copy(org);
    return *this;
}

template<class T>
inline CON_Matrix_t<T>& CON_Matrix_t<T> :: operator=(const T& org)
{
    if (exist())
    {
        unsigned int i,j;
        for (i=0;i<rows;i++)
            for (j=0;j<columns;j++)
                (* this)[i][j] = org;
    }
    return *this;
}

// define addition and substraction operators

template<class T>
inline CON_Matrix_t<T>& CON_Matrix_t<T> :: operator+=(const CON_Matrix_t<T>& org)
{
    static CON_Matrix_t<T> errorVar;
// test size of the source CON_Matrix_t
    if (org.rows!=rows || org.columns!=columns || !org.data || !data)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator+=","Dimension mismatch!",ut_error_array_dimmatch);
    }
// copy data
    int i,j;
    for (i=0;i<org.rows;i++)
        for (j=0;j<org.columns;j++)
            (*this)[i][j] += org[i][j];

    return *this;
}

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator+(const CON_Matrix_t<T>& org)
{
// create target object
    CON_Matrix_t<T> target;
// test size of the source CON_Matrix_t
    if (org.rows!=rows || org.columns!=columns || !org.data || !data)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator+","Dimension mismatch!",ut_error_array_dimmatch);
        return target;
    }

    try
    {
        target.allocate(rows,columns);
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator+","Allocation failed!",ut_error_array_alloc);
        return CON_Matrix_t<T>();
    }
// make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
// copy data
    int i,j;
    for (i=0;i<org.rows;i++)
        for (j=0;j<org.columns;j++)
            target[i][j] = (*this)[i][j] + org[i][j];

    return target;
}

template<class T>
inline CON_Matrix_t<T>& CON_Matrix_t<T> :: operator-=(const CON_Matrix_t<T>& org)
{
    static CON_Vector_t<T> errorVar;
// test size of the source CON_Matrix_t
    if (org.rows!=rows || org.columns!=columns || !org.data || !data)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator-=","Dimension mismatch!",ut_error_array_dimmatch);
    }
// copy data
    int i,j;
    for (i=0;i<org.rows;i++)
        for (j=0;j<org.columns;j++)
            (*this)[i][j] -= org[i][j];

    return *this;
}

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator-(const CON_Matrix_t<T>& org)
{
// create target object
    CON_Matrix_t<T> target;
// test size of the source CON_Matrix_t
    if (org.rows!=rows || org.columns!=columns || !org.data || !data)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator-","Dimension mismatch!",ut_error_array_dimmatch);
        return target;
    }
    try
    {
        target.allocate(rows,columns);
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator-","Allocation failed!",ut_error_array_alloc);
        return CON_Matrix_t<T>();
    }
// make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
// copy data
    int i,j;
    for (i=0;i<org.rows;i++)
        for (j=0;j<org.columns;j++)
            target[i][j] = (*this)[i][j] - org[i][j];

    return target;
}

// CON_Matrix_t multiplication

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator*(const CON_Matrix_t<T>& org) const
{
// Create target object
    CON_Matrix_t<T> target;

    multiply(*this, org, target);

    return target;
}

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator%(const CON_Matrix_t<T>& org) const
{
    if (org.height() != this->height() || org.length() != this->length())
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator%","Dimension mismatch!",ut_error_array_dimmatch);

    }
    CON_Matrix_t<T> target(org.height(),org.length());
    for (int i =0;i<target.height();i++) for (int j=0;j<target.length();j++)
        target[i][j] = org[i][j] * (*this)[i][j];

    return target;
}

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator/(const CON_Matrix_t<T>& org) const
{
    if (org.height() != this->height() || org.length() != this->length())
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator/","Dimension mismatch!",ut_error_array_dimmatch);

    }
    CON_Matrix_t<T> target(org.height(),org.length());
    try
    {
        for (int i =0;i<target.height();i++) for (int j=0;j<target.length();j++)
            target[i][j] = (*this)[i][j] / org[i][j];
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator/","Division by zero!",ut_error_unknown);

    }

    return target;
}


template<class T>
inline void CON_Matrix_t<T> :: multiply(const CON_Matrix_t<T>& in1, const CON_Matrix_t<T>& in2, CON_Matrix_t<T>& out)
{
// Test size of matrices
    if (in1.length() != in2.height())
    {
//        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiply(...)","Dimension mismatch!",ut_error_array_dimmatch);
        return;
    }

    try
    {
        out.allocate(in1.height(), in2.length());
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiply(...)","Allocation failed!",ut_error_array_alloc);
        return;
    }

// Copy data
    long i, j, k;
    T aux1 = 0;
    T aux2 = 0;
    const long nOutRows = out.height();
    const long nOutColumns = out.length();
    const long nIn1Columns = in1.length();
    for (i = 0; i < nOutRows; i++)
        for (j = 0; j < nOutColumns; j++)
        {
            out[i][j] = 0;
            for (k = 0; k < nIn1Columns; k++)
            {
                aux1 = in1[i][k];
                aux2 = in2[k][j];
                if (!aux1 || !aux2)  // if one operand is zero
                    continue;

                if (aux1 == 1) // if first operand is one
                {
                    out[i][j] += aux2;
                    continue;
                }

                if (aux2 == 1) // if second operand is one
                {
                    out[i][j] += aux1;
                    continue;
                }
                if (aux1 == -1) // if first operand is -one
                {
                    out[i][j] -= aux2;
                    continue;
                }

                if (aux2 == -1) // if second operand is -one
                {
                    out[i][j] -= aux1;
                    continue;
                }

                out[i][j] += aux1 * aux2;
            }
        }

    return;
}


// Multiplication of Matrix with diagonal Matrix on the right side, diagonal values are stored in a vector

template<class T>
inline void CON_Matrix_t<T> :: multiplyrightdiagonal(const CON_Vector_t<T>& in,CON_Matrix_t<T>& out)
{
    long i,j;

// Test size of matrices
    if ((*this).length() != in.length())
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiplyrightdiagonal(...)","Dimension mismatch!",ut_error_array_dimmatch);
        return;
    }

    try
    {
        out.allocate((*this).height(), in.length());
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiplyrightdiagonal(...)","Allocation failed!",ut_error_array_alloc);
        return;
    }

    for ( i = 0; i < out.height(); i ++)
        for ( j= 0; j< out.length(); j ++)
              out[i][j] = (*this)[i][j]* in[j];


}

// multiplications with a diagonalmatrix on the left side, diagonal matrix values in vector in1
template<class T>
inline void CON_Matrix_t<T> :: multiplyleftdiagonal(const CON_Vector_t<T>& in, CON_Matrix_t<T>& out)
{
    long i,j;

   // Test size of matrices
    if (in.length() != (*this).height())
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiplyleftdiagonal(...)","Dimension mismatch!",ut_error_array_dimmatch);
        return;
    }

    try
    {
        out.allocate(in.length(), (*this).length());
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiplyleftdiagonal(...)","Allocation failed!",ut_error_array_alloc);
        return;
    }


   for ( i = 0; i < out.height(); i ++)
       for ( j= 0; j< out.length(); j ++)
              out[i][j] = (*this)[i][j]* in[i];


}


// multiplication with CON_Vector_t

template<class T>
inline CON_Vector_t<T> CON_Matrix_t<T> :: operator*(const CON_Vector_t<T>& org) const
{
// create target object
    CON_Vector_t<T> target;
// test size of matrices
    if (org.length()!=columns)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator*","Wrong dimensions!",ut_error_array_dimmatch);
        return target;
    }
    try
    {
        target.allocate(rows);
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator*","Allocation failed!",ut_error_array_alloc);
        return CON_Vector_t<T>();
    }
// copy data
    int i,k;
    for (i=0;i<target.length();i++)
    {
        target[i] = 0;
        for (k=0;k<columns;k++)
            target[i] += (*this)[i][k] * org[k];
    }

    return target;
}


template<class T>
inline void CON_Matrix_t<T> :: multiply(const CON_Matrix_t<T>& in1, const CON_Vector_t<T>& in2, CON_Vector_t<T>& out)
{
// Test size of matrices
    if (in1.length() != in2.length())
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiply(...)","Wrong dimensions!",ut_error_array_dimmatch);
        return;
    }

    try
    {
        out.allocate(in1.height());
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::multiply(...)","Allocation failed!",ut_error_array_alloc);
        return;
    }

// Copy data
    long i, k;
    T aux1 = 0;
    T aux2 = 0;
    const long nOutRows = out.length();
    const long nIn1Columns = in1.length();
    for (i = 0; i < nOutRows; i++)
    {
        out[i] = 0;
        for (k = 0; k < nIn1Columns; k++)
        {
            aux1 = in1[i][k];
            aux2 = in2[k];
            if (!aux1 || !aux2)  // if one operand is zero
                continue;

            if (aux1 == 1) // if first operand is one
            {
                out[i] += aux2;
                continue;
            }

            if (aux2 == 1) // if second operand is one
            {
                out[i] += aux1;
                continue;
            }
            if (aux1 == -1) // if first operand is -one
            {
                out[i] -= aux2;
                continue;
            }

            if (aux2 == -1) // if second operand is -one
            {
                out[i] -= aux1;
                continue;
            }

            out[i] += aux1 * aux2;
        }
    }
    return;
}


// multiplication with scalar

template<class T>
inline CON_Matrix_t<T>& CON_Matrix_t<T> :: operator*=(const T& org)
{
    if (exist())
    {
        int i,j;
        for (i=0;i<rows;i++)
            for (j=0;j<columns;j++)
                (*this)[i][j] *= org;
    }

    return *this;
}

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator*(const T& org) const
{
// create target object
    CON_Matrix_t<T> target;
    CON_Matrix_t<T> error;
    try
    {
        target.allocate(rows,columns);
    }
    catch (...)
    {
        ASSERT(FALSE);
        throw Error("CON_Matrix_t::operator*","Allocation failed!",ut_error_array_alloc);
        error.allocate(0,0);
        return error;
    }
// make sure, the copy constructor leaves the allocated memory alone
    target.ownmem = false;
    target.reset_ownmem = true;
// copy data
    int i,j;
    for (i=0;i<rows;i++)
        for (j=0;j<columns;j++)
            target[i][j] = (*this)[i][j] * org;

    return target;
}

// equality

template<class T>
inline bool CON_Matrix_t<T> :: operator==(const CON_Matrix_t<T>& org) const
{
// if both matrices are empty, the ARE equal
    if (!exist() && !org.exist())
        return true;
// if the matrices have different dimensions, the are not equal
    if (rows != org.rows)
        return false;
    if (columns != org.columns)
        return false;
// check for equality of the elements
    int i,j;
    for(i=0;i<rows;i++)
        for (j=0;j<columns;j++)
            if ((*this)[i][j] != org[i][j])

                return false;
    return true;
}

template<class T>
inline bool CON_Matrix_t<T> :: operator!=(const CON_Matrix_t<T>& org) const
{
    return !((*this)==org);
}

template<class T>
inline CON_Vector_t<T> CON_Matrix_t<T> :: column(int index) const
{
    CON_Vector_t<T> target;
// if the matrix has no column of this index
    if (index < 0 || index >= columns)
        return target;

    target.allocate((*this).rows);
    for (int i=0;i<(*this).rows;i++)
        target[i] = (*this)[i][index];

    return target;
}

template<class T>
inline void CON_Matrix_t<T> :: setrow(const CON_Vector_t<T>& org, int index)
{
// if both matrices are empty
    if (!exist() || !org.exist())
        return;
// if the matrices have different dimensions
    if (columns != org.length())
        return;
// if the matrix has no row of this index
    if (index < 0 || index >= rows)
        return;
    for(int i=0;i<columns;i++)
        (*this)[index][i] = org[i];
}

template<class T>
inline void CON_Matrix_t<T> :: setcolumn(const CON_Vector_t<T>& org, int index)
{
// if both matrices are empty
    if (!exist() || !org.exist())
        return;
// if the matrices have different dimensions
    if (rows != org.length())
        return;
// if the matrix has no column of this index
    if (index < 0 || index >= columns)
        return;
    for(int i=0;i<rows;i++)
        (*this)[i][index] = org[i];
}

template<class T>
inline void CON_Matrix_t<T> :: removerow(int index)
{
    ASSERT(index >= 0 && index < rows);
    if (index < 0 || index >= rows)
        return;

    CON_Matrix_t<T> target(rows - 1, columns);
    for (int i = 0; i < rows; ++i)
        if (i != index)
            target.setrow((*this)[i], (i < index ? i : i - 1));

    *this = target;
}

template<class T>
inline void CON_Matrix_t<T> :: removerows(int index, int number)
{
    int i,c;
    for (i=index,c=0;c<number;i++,c++)
        removerow(i);
}

template<class T>
inline void CON_Matrix_t<T> :: removecolumn(int index)
{
    ASSERT(index >= 0 && index < columns);
    if (index < 0 || index >= columns)
        return;

    CON_Matrix_t<T> target(rows, columns - 1);
    for (int i = 0; i < columns; ++i)
        if (i != index)
            target.setcolumn(column(i), (i < index ? i : i - 1));

    *this = target;
}

template<class T>
inline void CON_Matrix_t<T> :: removecolumns(int index, int number)
{
    int i,c;
    for (i=index,c=0;c<number;i++,c++)
        removecolumn(i);
}

// frobenius norm

template<class T>
inline double CON_Matrix_t<T> :: norm() const
{
// if CON_Matrix_t is empty, return 0
    if (!exist())
        return 0.0;

// compute sum of squares of elements
    double N=0;
    int i,j;
    for (i=0;i<rows;i++)
        for (j=0;j<columns;j++)
            N += (double)(*this)[i][j] * (double)(*this)[i][j];

    return N;
}

// sum over all elements
template<class T>
inline double CON_Matrix_t<T> :: sum() const
{
// if CON_Matrix_t is empty, return 0
    if (!exist())
        return 0.0;

// compute sum of elements
    double N=0;
    int i,j;
    for (i=0;i<rows;i++)
        for (j=0;j<columns;j++)
            N += (double)(*this)[i][j];

    return N;
}

// mean over all elements
template<class T>
inline double CON_Matrix_t<T> :: mean() const
{
// if CON_Matrix_t is empty, return 0
    if (!exist())
        return 0.0;

    return sum()/(double)(rows*columns);
}

// Gauss

#ifdef WIN32
template<class T>
inline bool CON_Matrix_t<T> :: gauss(std::valarray<long double> & x, const std::valarray<long double> & _b) const
{
    CON_Matrix_t<T> A = *this;
    std::valarray<long double> b = _b;

    for (int row = 0; row < A.height()-1; row++)
    {
        int maxval = A.getMaxPosition(row);
        if (Abs(A[maxval][row]) < 1e-50)
        {
        // Singular matrix solving set of linear equations:
            return false;
        }

        if (maxval != row)
        {
            A.swaprows (row, maxval);
            swap (b[row], b[maxval]);
        }

        for (int i = row+1; i < A.height(); i++)
        {
            if (A[i][row] != 0)
            {
                double coeff = A[i][row]/A[row][row];
                A[i][row] = 0;
                for (int j = row+1; j < A.length(); j++)
                {
                    A[i][j] -= coeff * A[row][j];
                }
                b[i] -= coeff * b[row];
            }
        }
    }

    int N = A.length() - 1;

    x[N] = b[N]/A[N][N];
    for (int i = N-1; i >= 0; i--)
    {
        long double sum = 0;
        for (int j = i+1; j <= N; j++)
        {
            sum += A[i][j]*x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return true;      // Successful operation
}
#endif

// transpose CON_Matrix_t

template<class T>
void CON_Matrix_t<T> :: transpose()
{
    if (rows==columns) // in the symmetric case, we can transpose without actually copying the data
    {
        int i,j;
        T   a;
        for (i=0;i<rows;i++)
            for (j=0;j<i;j++)
            {
                a=(*this)[i][j];
                (*this)[i][j]=(*this)[j][i];
                (*this)[j][i]=a;
            }
    }
    else // for the unsymmetric case an intermediate copy of the data is necessary
    {
    // copy data, force copying of data CON_Block_t
        bool     temp1 = ownmem;
        unsigned temp2 = size;
        ownmem = true;
        CON_Matrix_t<T> aux1;
        aux1=(*this);
        ownmem = temp1;
    // create new object
    // change by Steinstraeter: If size != rows*columns, reusing data (in the way 'else' handles this)
    //                          destroys the data block (typically owned by another instance, see assign_reduced()).
        if (ownmem || (size != rows*columns))
        {
            clear();
            try
            {
                allocate(aux1.columns,aux1.rows);
            }
            catch (...)
            {
                ASSERT(FALSE);
                throw Error("CON_Matrix_t::transpose()","Allocation failed!",ut_error_array_alloc);

                return;
            }
        }
        else
        {
            T *aux2 = data;
            clear();
            assign(aux1.columns,aux1.rows,aux2);
            size = temp2;
        }

    // copy data
        int i,j;
        for (i=0;i<rows;i++)
            for (j=0;j<columns;j++)
                (*this)[i][j] = aux1[j][i];
    }
}

// mirror CON_Matrix_t

template<class T>
inline void CON_Matrix_t<T> :: mirror_horizontal()
{
// if CON_Matrix_t is empty, do nothing
    if (exist())
    {
        int i,j;
        T   aux;
        for (i=0;i<rows;i++)
        {
            for (j=0;j<columns-j-1;j++)
            {
                aux = (*this)[i][j];
                (*this)[i][j]=(*this)[i][columns-j-1];
                (*this)[i][columns-j-1]=aux;
            }
        }
    }
}

template<class T>
inline void CON_Matrix_t<T> :: mirror_vertical()
{
// if CON_Matrix_t is empty, do nothing
    if (exist())
    {
        int i,j;
        T   aux;
        for (i=0;i<rows-i-1;i++)
        {
            for (j=0;j<columns;j++)
            {
                aux = (*this)[i][j];
                (*this)[i][j]=(*this)[rows-i-1][j];
                (*this)[rows-i-1][j]=aux;
            }
        }
    }
}

template<class T>
inline void CON_Matrix_t<T> :: swaprows(int row1, int row2)
{
    for (int col = 0; col < columns; col++)
    {
        swap((*this)[row1][col], (*this)[row2][col]);
    }
}

// concatination of two matrices (horizontal, join columns)

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator|(const CON_Matrix_t<T>& org) const
{
    CON_Matrix_t<T> target;
    if (exist() || org.exist())
    {
    // test size of matrices
        if (org.rows!=rows || !org.data || !data)
        {
            ASSERT(FALSE);
            throw Error("CON_Matrix_t::operator|","Wrong dimensions!",ut_error_array_dimmatch);
            return target;
        }
    // create target object
        try
        {
            target.allocate(rows,columns+org.columns);
        }
        catch (...)
        {
            ASSERT(FALSE);
            throw Error("CON_Matrix_t::operator|","Allocation failed!",ut_error_array_alloc);
            return CON_Matrix_t<T>();
        }
    // make sure, the copy constructor leaves the allocated memory alone
        target.ownmem = false;
        target.reset_ownmem = true;
    // fill in values
        int i,j,k;
        for (i=0;i<rows;i++)
        {
//            for (j=k=0;j<columns;j++,k++)
            for (j=0,k=0;j<columns;j++,k++)
                target[i][k] = (*this)[i][j];
            for (j=0;j<org.columns;j++,k++)
                target[i][k] = org[i][j];
        }
    }

    return target;
}

// concatination of two matrices (vertical, join rows)

template<class T>
inline CON_Matrix_t<T> CON_Matrix_t<T> :: operator&(const CON_Matrix_t<T>& org)  const
{
    CON_Matrix_t<T> target;
    if (exist() || org.exist())
    {
    // test size of matrices
        if (org.columns!=columns || !org.data || !data)
        {
            ASSERT(FALSE);
            throw Error("CON_Matrix_t::operator&","Wrong dimensions!",ut_error_array_dimmatch);
            return target;
        }
        // create target object
        try
        {
            target.allocate(rows+org.rows,columns);
        }
        catch (...)
        {
            ASSERT(FALSE);
            throw Error("CON_Matrix_t::operator&","Allocation failed!",ut_error_array_alloc);
            return CON_Matrix_t<T>();
        }
    // make sure, the copy constructor leaves the allocated memory alone
        target.ownmem = false;
        target.reset_ownmem = true;
    // fill in values
        int i,j,k;
        for (j=0;j<columns;j++)
        {
//            for (i=k=0;i<rows;i++,k++)
            for (i=0,k=0;i<rows;i++,k++)
                target[k][j] = (*this)[i][j];
            for (i=0;i<org.rows;i++,k++)
                target[k][j] = org[i][j];
        }
    }
    return target;
}

// harakiri function: do not use unless authorised
template<class T>
inline T* CON_Matrix_t<T> :: getdata()
{
    return data;
}

template<class T>
inline const T* CON_Matrix_t<T> :: getdata() const
{
    return data;
}

// added by Steinstraeter: 'f' if T = float, 'd' if T = double, '?' otherwise.
template<class T>
inline char CON_Matrix_t<T>::getElementType()
{
  return CON_Vector_t<T>::getElementType();
}

// added by Steinstraeter:
//   string conversion, matlab notation
template<class T>
std::string CON_Matrix_t<T>::toString(int maxLen, int setw_value, char format, int setprecision_value, bool use_cr) const
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
      if (maxLen < 0)
        maxLen = height()*length();

      bool abort = false;
      int count = 0;
      for (int i = 0; i < height(); i++)
        {
          for (int j = 0; j < length()-1; j++)
            {
              if (count >= maxLen)
                {
                  abort = true;
                  break;
                }
              s << std::setw(setw_value) << (*this)[i][j] << " ";
              count++;
            }
          if (abort)
            break;
          if (count >= maxLen)
            {
              abort = true;
              break;
            }
          s << std::setw(setw_value) << (*this)[i][length()-1];
          count++;
          if (i < height()-1) {
            if (use_cr)
              s << "; ...\n ";
            else
              s << "; ";
          }
        }
      if (abort)
        {
          int furtherElements = height()*length()-count;
          char pl[] = "s";
          if (furtherElements == 1)
            pl[0] = 0;
          s << "... <" << furtherElements << " further element" << pl <<", "<< height() << "x" << length() << " matrix>";
        }
    }
  s << "]";

  return s.str();
}
// added by Steinstraeter: save matrix elements to binary file
template<class T>
inline bool CON_Matrix_t<T>::save(const char *filename) const
{
  char outputFormat = 'd';
  if (getElementType() == 'f')
    outputFormat = 'f';
  std::ofstream *fptr = CON_Vector_t<T>::save_header(filename, outputFormat, height(), length());

  if (*fptr)
    for (int i = 0; i < height() && (*fptr); i++)
      (*this)[i].rawsave(*fptr, outputFormat);

  bool success = (*fptr? true:false);
  fptr->close();
  delete fptr;
  return success;
}

// added by Steinstraeter: read vector elements from binary file
template<class T>
inline bool CON_Matrix_t<T>::load(const char *filename)
{
  char fileFormat;
  int fileRows, fileColumns;
  std::ifstream *fptr = CON_Vector_t<T>::read_header(filename, fileFormat, fileRows, fileColumns);

  if (fileFormat == 'd' || fileFormat == 'f')
    {
      allocate_or_clear(fileRows, fileColumns);
      for (int i = 0; i < height() && (*fptr); i++)
        (*this)[i].rawread(*fptr, fileFormat);
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
inline void CON_Matrix_t<T> :: init()
{
    data = NULL;
    datarows = NULL;
    rows = columns = size = 0;
    ownmem = false;
    reset_ownmem = false;
}

// common code used in assignment operator and copy constructor

template<class T>
void CON_Matrix_t<T> :: copy(const CON_Matrix_t<T>& org)
{
// first empty the object and free all memory
    clear();
// copy the dimensions and memory flags
    rows = org.rows;
    columns = org.columns;
    if (org.reset_ownmem)
        ownmem = true;
    else
        ownmem = org.ownmem;

  // removed by Steinstraeter: org.size != org.rows * org.columns is possible
    //size = org.size;

  if (org.ownmem || !org.reset_ownmem) // TRK11022002, if not ownmem=false AND reset_ownmem=true,
                                           // as used in temporary objects that are put on stack by operator functions

    {
    // added by Steinstraeter: org.size != org.rows*org.columns is possible
      size = rows*columns;

    // try to allocate memory, if unsuccessful, throw exception
        if (org.size)
            try
            {
        // changed by Steinstraeter: org.size != org.rows*org.columns is possible
                //   data = new T[org.size];
                data = new T[size];
            }
            catch (std::bad_alloc&)
            {
                ASSERT(FALSE);
                throw Error("CON_Matrix_t::copy()","Allocation failed!",ut_error_array_alloc);
                return;
            }

        if (org.rows)
            try
            {
                datarows = new CON_Vector_t<T>[org.rows];
            }
            catch (std::bad_alloc&)
            {
                ASSERT(FALSE);
                throw Error("CON_Matrix_t::copy()","Allocation failed!",ut_error_array_alloc);
                return;
            }
    // copy data
        unsigned int ii;
    // changed by Steinstraeter: BGN
    //  org.size != org.rows*org.columns is possible
        //      for (ii=0;ii<org.rows;ii++)
        //          datarows[ii].assign(org.columns,data+org.hashtable[ii]);
      //        memcpy(data,org.data,size*sizeof(T));  //TRK11022002
    hashtable.allocate(rows);
        for (ii=0;ii<org.rows;ii++)
      {
          hashtable[ii] = ii * columns;
          memcpy(data+hashtable[ii], org.data+org.hashtable[ii], columns*sizeof(T));
              datarows[ii].assign(columns,data+hashtable[ii]);
      }
    // changed by Steinstraeter: END
        ownmem = true;                           //TRK11022002
        // for (ii=0;ii<org.size;ii++)
        //     data[ii] = org.data[ii];
    }
    else // if the original does not own its data
       // comment added by Steinstraeter: Even if org does not own its data the data are normally copied.
       //                                 The data are not copied iff org.ownmem == false and org.reset_ownmem == true.
       //                                 org indicates with these parameters, that the responsibility for the data should be
       //                                 transfered to *this.
    {
    // added by Steinstraeter:
      size = org.size;

        data = org.data;
    // try to allocate memory, if unsuccessful, throw exception
        if (org.rows)
            try
            {
                datarows = new CON_Vector_t<T>[org.rows];
            }
            catch (std::bad_alloc&) // seems to be Borland specific (otherwise bad_alloc)
            {
                ASSERT(FALSE);
                throw Error("CON_Matrix_t::copy()","Allocation failed!",ut_error_array_alloc);
                return;
            }
    // copy data
        for (unsigned int ii=0;ii<org.rows;ii++)
            datarows[ii].assign(org.columns,data+org.hashtable[ii]);

      hashtable=org.hashtable;
    }

  // removed by Steinstraeter: May be incorrect if org.size != org.rows * org.columns
    //hashtable=org.hashtable;
}

template <class T>
T CON_Matrix_t<T>::getMaxi() const
{
    T tMaxi = 0;
    if (exist())
    {
        tMaxi = datarows[0].getMaxi();
        for (long nSample = 1; nSample < rows; ++nSample)
            tMaxi = max(tMaxi, datarows[nSample].getMaxi());
    }

    return tMaxi;
}

template <class T>
T CON_Matrix_t<T>::getMini() const
{
    T tMini = 0;
    if (exist())
    {
        tMini = datarows[0].getMini();
        for (long nSample = 1; nSample < rows; ++nSample)
            tMini = min(tMini, datarows[nSample].getMini());
    }

    return tMini;
}
/////////////////////////////////////////////////////////////////////////////
// CON_Matrix_t Serialization
//
#ifdef WIN32
template<class T>
inline void CON_Matrix_t<T>::Serialize(class CON_AbstractPersistenceArchive_i& ar)
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

        ar << rows;
        ar << columns;
        ar << size;
        ASSERT(size == rows * columns);
        ar.Write(data, columns * rows  * sizeof(T));
    }
    else
    {
        ar >> rows;
        ar >> columns;
        ar >> size;
        ASSERT(size == rows * columns);

        allocate(rows, columns);
        reset_ownmem = false;

        ar.Read(data, rows * columns * sizeof(T));
    }
}

template<class T>
inline CON_AbstractPersistenceArchive_i& operator<<(CON_AbstractPersistenceArchive_i& ar, const CON_Matrix_t<T>& matrix)
{
    const_cast< CON_Matrix_t<T>& >(matrix).Serialize(ar);

    return ar;
}

template<class T>
inline CON_AbstractPersistenceArchive_i& operator>>(CON_AbstractPersistenceArchive_i& ar, CON_Matrix_t<T>& matrix)
{
    matrix.Serialize(ar);

    return ar;
}
#endif
#endif // __utMatrix_t_H__
