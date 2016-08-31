//$7    17.09.2002    TRK added string interface
//$6    12.10.2000  Matthias D. Removed capital letters from file names
//$5    09.08.1999  TRK     use namespace std for BC++ as well
//$4    10.05.1999  Frank N. Changed includes
//$3    05.05.1999  Frank N. Changed include concept
//$2    03.05.1999    Frank N. submitted for TRK
//$1    12.04.1999    Frank N.    split from ut_array_t.cpp
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
/*-------------------------------------------------------------------------*/
/* File name : Error.cpp                                                  */
/*-------------------------------------------------------------------------*/

#include "StdUtilitiesPch.h"
#include <string>            // STL string class
// typedef std::basic_string<char> string;

using namespace std;

#include "Error.h"

// Exception class
//------------------------------------------------------------------------------

// constructor: create exception with description of origin and content


Error :: Error(const char* origin,const char* description,const unsigned c)
 {
  // check, if error code is valid
  if (c <= ut_error_mincode || c >= ut_error_maxcode) code = ut_error_unknown;
  else                                                code = c;
  // allocate strings (if error occurs, throw on)
  try {message  = description;} catch(...) {throw;}
  try {location = origin;}      catch(...) {throw;}
 }

Error :: Error(const std::string &origin,const std::string &description,const unsigned c)
 {
  // check, if error code is valid
  if (c <= ut_error_mincode || c >= ut_error_maxcode) code = ut_error_unknown;
  else                                                code = c;
  // allocate strings (if error occurs, throw on)
  try {message  = description;} catch(...) {throw;}
  try {location = origin;}      catch(...) {throw;}
 }

// constructor: copy exception

Error :: Error(const Error& org)
 {
  code = org.code;
  // allocate strings (if error occurs, throw on)
  try {message  = org.message;}  catch(...) {throw;}
  try {location = org.location;} catch(...) {throw;}
 }

// destructor: destroy exception object

Error :: ~Error()
 {
 }

// return information "where" an exception occurred

const char * Error :: where()
 {
  return location.c_str();
 }

const string Error :: Where()
 {
  return location;
 }

// return information "what" has happened

const char * Error :: what()
 {
  return message.c_str();
 }

const string Error :: What()
 {
  return message;
 }
