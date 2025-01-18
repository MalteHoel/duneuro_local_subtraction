// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_GMSH_TENSOR_READER_HH
#define DUNEURO_GMSH_TENSOR_READER_HH

#include <fstream>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/timer.hh>

#include <duneuro/io/data_tree.hh>

namespace duneuro
{
  template <class G>
  class GmshTensorReader
  {
  public:
    typedef typename G::ctype ctype;
    enum { dim = G::dimension };
    typedef Dune::FieldMatrix<ctype, dim, dim> TensorType;

    static void read(const std::string& filename, std::vector<TensorType>& tensors,
                     DataTree dataTree = DataTree())
    {
      Dune::Timer timer;
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "tensor file " << filename << " could not be opened");
      }
      // assuming isotropic conductivities
      double value;
      while (stream >> value) {
        tensors.push_back(isotropicTensor(value));
      }
      if (!stream.eof()) {
        DUNE_THROW(Dune::IOError, "not all tensors could be read. check the file formatting");
      }
      dataTree.set("tensors", tensors.size());
      dataTree.set("time", timer.elapsed());
    }

  private:
    static TensorType isotropicTensor(ctype value)
    {
      TensorType t;
      for (unsigned int r = 0; r < t.N(); ++r) {
        for (unsigned int c = 0; c < t.M(); ++c) {
          t[r][c] = r == c ? value : 0.0;
        }
      }
      return t;
    }
  };
}

#endif // DUNEURO_GMSH_TENSOR_READER_HH
