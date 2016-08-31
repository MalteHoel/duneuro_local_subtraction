//$24    21.05.2003  Anwander A. changed function PseudoInverse to PseudoInverseAndSolve
//$23    16.10.2002  Frank        Added tutest()
//$22    30.07.2002  Matthias D. Added PseudoInverse() and TSVD using float values
//$21    18.02.2002    Anwander A. round (double) defined in /usr/include/bits/mathcalls.h
//$20   11.10.2000  Matthias D. Removed capital letters from file names, disabled methods using valarray using ifdef WIN32
//$19    14.08.2000    Frank N.    For SVD() and PseudoInverse() use non const matrix since these functions change the matrix
//$18    20.03.2000    Frank N.    Added expj() and evaluate()
//$17    14.03.2000    Frank Z.    Added SolidAngle()
//$16    13.03.2000    Frank N.    Fixed const problem with Crossprod()
//$15    02.12.1999    Frank Z.    added round()
//$14    04.11.1999    Frank Z.    moved file to utilities
//$13    22.10.1999    Frank N.    Changed class names
//$12    21.10.1999    Frank Z.    moved functions to SignalSources
//$11    20.10.1999    Frank Z.    changed concept of SignalSources, major code revision
//$10    15.09.1999    Frank Z.    changed MatchElectodes
//$9    26.08.1999    Frank Z.    added pseudoinverse
//$8    17.08.1999    Frank Z.    added math. template functions
//$7    09.08.1999    TRK            Added #include <ut_datawrapper.h>
//$6    05.08.1999    Frank N.    Added forward declarations
//$5    02.08.1999    Frank N.    Moved to dataCore
//$4    22.07.1999    Frank N.    Second argument of MatchElectrodes() became const
//$3    15.07.1999    Frank Z.    added MaktchElectrodes() and Sort()
//$2    29.06.1999    Frank Z.    added MakeTriangulation()
//$1    11.06.1999    Frank Z.    created
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
// File name : ut-auxil.h
//-------------------------------------------------------------------------

#ifndef __UT_AUXIL_H__
#define __UT_AUXIL_H__

#include <complex>

#include "CON_Vector_t.h"
#include "CON_Matrix_t.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef pow10
#define pow10(arg) pow(10.0,arg)
#endif

// Forward declarations
#ifdef WIN32
template <class T>
class std::valarray;
#endif

template <class DataType>
void Cart2Sphere(const CON_Vector_t< DataType >& Cart,CON_Vector_t< DataType >& Sphere);

UTILITIES_EXPORT int     ShadingFactor(CON_Vector_t<float>& v,CON_Vector_t<float>& l,CON_Vector_t<float>& n);
UTILITIES_EXPORT float     Norm(CON_Vector_t<float>&);

template <class DataType>
DataType Scalar(const CON_Vector_t<DataType>& aa, const CON_Vector_t<DataType>& ab);

template <class DataType>
void Crossprod(CON_Vector_t<DataType>&, const CON_Vector_t<DataType>&, const CON_Vector_t<DataType>&);

UTILITIES_EXPORT long     round(float);
#ifdef WIN32
UTILITIES_EXPORT long     round(double);
#endif
UTILITIES_EXPORT void     shell(int n,CON_Vector_t<float>& a,CON_Vector_t<int>& index);    // sorting algorithm
UTILITIES_EXPORT int    Ludcmp(CON_Matrix_t<double>& a,int n,CON_Vector_t<int>& indx,int d);
UTILITIES_EXPORT void    Lubksb(CON_Matrix_t<double>& a,int n,CON_Vector_t<int>& indx,CON_Vector_t<double>& bb);
UTILITIES_EXPORT void    gaussj(CON_Matrix_t<double>& a, int n,CON_Matrix_t<double>& b, int m);
UTILITIES_EXPORT float     CrossSectionLinePlane(CON_Vector_t<float>&,CON_Vector_t<float>&,const CON_Vector_t<float>&,const CON_Vector_t<float>&);
UTILITIES_EXPORT int    TriEdgesInPlane(CON_Vector_t<float>& r0,CON_Vector_t<float>& r1,CON_Vector_t<float>& r2,CON_Vector_t<float>& n,CON_Vector_t<float>& p0,CON_Vector_t<float>& res);
UTILITIES_EXPORT int    SetToNewOrigin(CON_Matrix_t<float>& Pnts,int NumberPnts,CON_Vector_t<float>& Origin,CON_Vector_t<float>& NewX,CON_Vector_t<float>& NewY);
UTILITIES_EXPORT int    SetToOldOrigin(CON_Matrix_t<float>& Pnts,int NumberPnts,CON_Vector_t<float>& Origin,CON_Vector_t<float>& NewX,CON_Vector_t<float>& NewY);
UTILITIES_EXPORT void TriangleNorm(CON_Vector_t<float>& res,CON_Vector_t<float>& v1,CON_Vector_t<float>& v2,CON_Vector_t<float>& v3);
UTILITIES_EXPORT bool PseudoInverseAndSolve(int m,int n,int k, CON_Matrix_t<double>& a,const CON_Matrix_t<double>& b,CON_Matrix_t<double>& res);
UTILITIES_EXPORT bool PseudoInverse(CON_Matrix_t<double>& inMatrix, CON_Matrix_t<double>& outMatrix,double cond = 1E-6);
UTILITIES_EXPORT bool PseudoInverse(CON_Matrix_t<float>& inMatrix, CON_Matrix_t<float>& outMatrix,float cond = 1E-6);
UTILITIES_EXPORT int SVD(CON_Matrix_t<double>& a, CON_Vector_t<double>& w, CON_Matrix_t<double>& v);
UTILITIES_EXPORT int SVD(CON_Matrix_t<float>& a, CON_Vector_t<float>& w, CON_Matrix_t<float>& v);
UTILITIES_EXPORT int HasCommonEdge(CON_Matrix_t<unsigned int>& poly,int a,int b);
UTILITIES_EXPORT int CheckDiagonal(CON_Matrix_t<unsigned int>& poly,CON_Matrix_t<float>& points,int x,int y);
UTILITIES_EXPORT bool tutest(float ave1, float var1, unsigned long n1, float ave2, float var2, unsigned long n2, float *t, float *prob);
UTILITIES_EXPORT bool PseudoInverse_lapack(CON_Matrix_t<double>& inMatrix, CON_Matrix_t<double>& outMatrix);
UTILITIES_EXPORT bool Matrixinversion_lapack(CON_Matrix_t<double>& inMatrix, CON_Matrix_t<double>& outMatrix);
UTILITIES_EXPORT bool SVD_lapack(CON_Matrix_t<double>& inMatrix,CON_Matrix_t<double>& U,CON_Matrix_t<double>& D,CON_Matrix_t<double>& VT);
UTILITIES_EXPORT CON_Matrix_t<double> Multiply_Matrix_lapack(CON_Matrix_t<double> &A,CON_Matrix_t<double> &B);
UTILITIES_EXPORT CON_Vector_t<double> Multiply_Vector_Matrix(CON_Matrix_t<double>& inMatrix, CON_Vector_t<double>& vec);
// added by Steinstraeter: BGN
UTILITIES_EXPORT int svd_lapack(CON_Matrix_t<double>& A, CON_Matrix_t<double>& U, CON_Vector_t<double>& d, CON_Matrix_t<double>& Vt, CON_Vector_t<double> *w = NULL);
UTILITIES_EXPORT int svd_lapack(CON_Matrix_t<double>& A, CON_Vector_t<double>& d, CON_Vector_t<double> *w = NULL);
// added by Steinstraeter: END
/////////////////////////////////////////////////////////////////////////////
// Cart2Sphere()

//added by Dannhauer: begin
bool PseudoInverse_lapack2(CON_Matrix_t<double>& inMatrix, CON_Matrix_t<double>& outMatrix);
bool schur_lapack(CON_Matrix_t<double>& inMatrix, CON_Matrix_t<double>& outMatrix);
double epsi(double number);
//added by Dannhauer: end


template <class DataType>
void Cart2Sphere(const CON_Vector_t< DataType >& Cart, CON_Vector_t< DataType >& Sphere)
{
    DataType r, theta, alpha;

    r = sqrx(Cart[0] * Cart[0] + Cart[1] * Cart[1] + Cart[2] * Cart[2]);
    if (r > 1.e-19)
    {
        theta = acos(Cart[2] / r);
        if (fabs(Cart[0]) < 1.e-19)
        {
            if (Cart[1] > 0)
                alpha = M_PI / 2;
            else
                alpha = 3 * M_PI / 2;
            if (fabs(Cart[1]) < 1e-19)
                alpha = 0;
        }
        else
            alpha = atan2(Cart[1],Cart[0]);
    }
    else
    {
        theta = 0.0;
        alpha = 0.0;
    }
    Sphere[0] = r;
    Sphere[1] = theta;
    Sphere[2] = alpha;
}

/////////////////////////////////////////////////////////////////////////////
// Scalar()

template <class DataType>
inline DataType Scalar(const CON_Vector_t<DataType>& aa, const CON_Vector_t<DataType>& ab)
/******************************************************************************
dot product of two vectors aa and ab.
******************************************************************************/

{
  return aa * ab;
}

/////////////////////////////////////////////////////////////////////////////
// Crossprod()

template <class DataType>
inline void Crossprod(CON_Vector_t<DataType>& V3Out,const CON_Vector_t<DataType>& V1In, const CON_Vector_t<DataType>& V2In)
/******************************************************************************
CON_Vector_t product of two vectors aa and ab. V3 = V1 x V2
******************************************************************************/
{
     V3Out = V1In / V2In;
}

template <class DataType>
inline DataType SolidAngle(const CON_Vector_t<DataType>& origin,const CON_Vector_t<DataType>& x1,const CON_Vector_t<DataType>& x2,const CON_Vector_t<DataType>& x3)
// computation of the normalized solid angle, seen from the origin to
// triangle(r1,r2,r3).
{
    CON_Vector_t<DataType>  r1(3,0),r2(3,0),r3(3,0),r4(3,0);
    DataType n1, n2, n3, n4, n5, angle;

    for (int i=0;i<3;i++)
    {
        r1[i] = x1[i] - origin[i];
        r2[i] = x2[i] - origin[i];
        r3[i] = x3[i] - origin[i];
    }
    r4[0] = r2[1] * r3[2] - r2[2] * r3[1]; /* r2 cross r3 */
    r4[1] = r2[2] * r3[0] - r2[0] * r3[2];
    r4[2] = r2[0] * r3[1] - r2[1] * r3[0];
    n1 = r1.norm();
    n2 = r2.norm();
    n3 = r3.norm();
    n4 = Scalar(r1,r4);                   /* r1 * r4 */
    n5 = n1 * n2 * n3 + Scalar(r1,r2)*n3 + Scalar(r1,r3)*n2 + Scalar(r2,r3)*n1;
    angle = 1E-30;
    if ( fabs(n5) > 1E-30 )
        angle = 2 * atan2(n4, n5);
    if (angle<0)
        angle = 4*M_PI+angle;

    return angle;
}

template <class T>
inline T Abs (const T & x)
{
    return x > 0 ? x : -x;
}

template <class T>
inline std::complex<T> expj(T theta)
{
    return std::complex<T>(cos(theta), sin(theta));
}

template <class T>
inline T hypot(std::complex<T> z)
{
    return ::hypot(z.imag(), z.real());
}

template <class T>
inline T atan2(std::complex<T> z)
{
    return ::atan2(z.imag(), z.real());
}
#ifdef WIN32
template <class T>
inline std::complex<T> evaluate(const std::valarray< std::complex<T> >& topco, int nz, const std::valarray< std::complex<T> >& botco, int np, std::complex<T> z)
{
// evaluate response, substituting for z:
    return eval(topco, nz, z) / eval(botco, np, z);
}

template <class T>
inline std::complex<T> eval(const std::valarray< std::complex<T> >& coeffs, int npz, std::complex<T> z)
{
// evaluate polynomial in z, substituting for z:
    std::complex<T> sum = std::complex<T>(0.0);
    for (int i = npz; i >= 0; i--)
        sum = (sum * z) + coeffs[i];

    return sum;
}
#endif
template <class T>
inline std::complex<T> csqrt(std::complex<T> x)
{
    T r = hypot(x);
    std::complex<T> z = std::complex<T>(sqrx(0.5 * (r + x.real())),
    sqrx(0.5 * (r - x.real())));
    if (x.imag() < 0.0)
        z = std::complex<T>(z.real(), -z.imag());

    return z;
}

template <class T>
inline std::complex<T> sqr(std::complex<T> x)
{
    return x * x;
}

#endif /* ifndef __UT_AUXIL_H__ */
