//$10    17.09.2002    TRK added string interface
//$9    06.02.2001  Maurice     added ut_error_tree_paste
//$8    17.02.2000    TRK            added some error codes for BEM
//$7    09.08.1999  Frank N.    Use std::string
//$6    05.08.1999  TRK include cstring rausgeschmissen
//$5    10.05.1999  Frank N. Changed includes
//$4    05.05.1999  Frank N. Changed include concept
//$3    04.05.1999  FZ made functions exportable
//$2    03.05.1999  TRK            minor changes
//$1    12.04.1999    Frank N.    split from ut_array_t.h
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
// File name : Error.h
//-------------------------------------------------------------------------

#ifndef _ERROR
#define _ERROR

#include <string>            // STL string class

// Usually, general exceptions from the C++ exception hierarchy are caught at
// an early stage and thrown on as ASA exception. The infomation on the nature
// (what) and the location (where) of the exception should be concretised on
// each level of exception handling, until the problem can be finally treated.
//------------------------------------------------------------------------------

enum ut_errorcode{ut_error_mincode,                 // minimum possible code
                  //------------------------------------------------------------
                  // errors with array objects
                  ut_error_array_range,             // attempted access to an index out of range
                                                    // e.g. {vector<bool> a(3);bool i=a[55]; -> exception, index out of range! }
                  ut_error_array_alloc,             // memory could not be allocated
                  ut_error_array_index,             // attempt to allocate data with invalid dimensions (e.g. negative)
                  ut_error_array_dimmatch,          // dimensions of two arrays to be added, subtracted, or multiplied do not match as required for this operation
                                                    // e.g. {matrix<int> a(2,3),b(4,5);matrix<int> c = a*b; -> exception, columns of a and rows of b must agree! }
                  ut_error_array_selectrows,        // selection vector has a different length than full matrix
                                                    // in assign_reduced()
                  ut_error_array_selectcols,        // offset nd number of samples for column selection are not compatible
                                                    // in assign_reduced() or copy_reduced()
                  // errors with datawrapper
                  ut_error_geo_copy,                // object could not be copied (e.g. allocation fault or corrupted original)
                  ut_error_geo_alloc,
                  ut_error_geo_join,
                  ut_error_geo_dimension,           // attempted allocation with invalid dimensions (i.e. < 0)
                  ut_error_geo_timescale,           // attempted to define invalid time scale (e.g. end time < starting time)
                  ut_error_geo_fileio,              // file input/output error, e.g. failure of fclose() or fopen()
                  ut_error_geo_wronglength,         // vector has unexpected length
                  ut_error_geo_range,               // attempted access to an index out of range
                  ut_error_img_negspace,            // attempted to feed negative voxel spacing
                  ut_error_complist_type,            // attempted to access wrong compartmentlist type i.e. eeg_compute_compartment meg_compute_compartment render_compartment.
                  // errors with simulator
                  ut_error_wrongmodel,              // head model incorrect

                  // errors regarding the history tree
                  ut_error_tree_paste,                // could not paste feature
                  // errors for BSSPM
                  ut_error_bsspm,
                  ut_error_grid,
                  ut_error_gradiometer,
                  //------------------------------------------------------------
                  ut_error_maxcode,                 // maximum possible code
                  ut_error_unknown                  // unknown error

                };

class UTILITIES_EXPORT Error
{
public:
   Error(const char*,const char*,const unsigned);    // create exception
   Error(const std::string &,const std::string &, const unsigned);
   Error(const Error&);          // copy exception
   ~Error();                        // destroy exception
   const char *what();                 // return message string
   const char *where();                // return location string
   const std::string What();
   const std::string Where();
   unsigned code;                      // unique exception code
protected:
   std::string message;
   std::string location;
};

#endif    // _UT_ERROR
