// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_TEST_TEST_UTILITIES_HH
#define DUNEURO_TEST_TEST_UTILITIES_HH

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <fstream>

#include <dune/common/exceptions.hh>

namespace duneuro {

template<class T>
T norm(const std::vector<T>& vector) 
{
  return std::sqrt(std::inner_product(vector.begin(), vector.end(), vector.begin(), T(0.0)));
}

template<class T>
T relativeError(const std::vector<T>& test, const std::vector<T>& ref) 
{
  if(test.size() != ref.size()) {
    DUNE_THROW(Dune::RangeError, "size of test (" << test.size() << ") does not match size of ref (" << ref.size() << ")");
  }
  
  std::vector<T> diff;
  std::transform(test.begin(), test.end(), ref.begin(), std::back_inserter(diff),
                 [] (const T& test_val, const T& ref_val) {return (test_val - ref_val);});
  return norm(diff) / norm(ref);
}
  
template<class T>
T rdm(const std::vector<T>& test, const std::vector<T>& ref) 
{
  if(test.size() != ref.size()) {
    DUNE_THROW(Dune::RangeError, "size of test (" << test.size() << ") does not match size of ref (" << ref.size() << ")");
  }

  T norm_test = norm(test);
  T norm_ref = norm(ref);
  std::vector<T> diff;
  std::transform(test.begin(), test.end(), ref.begin(), std::back_inserter(diff),
                  [norm_test, norm_ref] (const T& test_val, const T& ref_val) 
                    {return (test_val/norm_test - ref_val/norm_ref);});
  return norm(diff);
}

template<class T>
T mag(const std::vector<T>& test, const std::vector<T>& ref) 
{
  if(test.size() != ref.size()) {
    DUNE_THROW(Dune::RangeError, "size of test (" << test.size() << ") does not match size of ref (" << ref.size() << ")");
  }

  return norm(test) / norm(ref);
}

template<class T>
void subtractMean(std::vector<T>& vector) 
{
  T mean = std::accumulate(vector.begin(), vector.end(), T(0.0)) / vector.size();
  for(T& entry : vector) {
    entry -= mean;
  }
  return;
}

template<class T>
std::vector<T> readFromTxt(const std::string& filename)
{
  std::ifstream inputStream(filename);
  if(!inputStream) {
    DUNE_THROW(Dune::IOError, "could not open " << filename << "!");
  }
  std::vector<T> valuesFromTxt;
  T value;
  while(inputStream >> value) {
    valuesFromTxt.push_back(value);
  }
  if(!inputStream.eof()) {
    DUNE_THROW(Dune::IOError, "reading file failed");
  }
  return valuesFromTxt;
}

} // namespace duneuro

#endif // DUNEURO_TEST_TEST_UTILITIES_HH
