#ifndef DUNEURO_GMSH_TENSOR_READER_HH
#define DUNEURO_GMSH_TENSOR_READER_HH

#include <fstream>

#include <dune/common/exceptions.hh>

#include <dune/common/fmatrix.hh>

namespace duneuro
{
  template <class G>
  class GmshTensorReader
  {
  public:
    typedef typename G::ctype ctype;
    enum { dim = G::dimension };
    typedef Dune::FieldMatrix<ctype, dim, dim> TensorType;

    static void read(const std::string& filename, std::vector<TensorType>& tensors)
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "tensor file " << filename << " could not be opened");
      }
      // assuming isotropic conductivities
      double value;
      while (stream >> value) {
        TensorType t;
        for (unsigned int r = 0; r < t.N(); ++r) {
          for (unsigned int c = 0; c < t.M(); ++c) {
            t[r][c] = r == c ? value : 0.0;
          }
        }
        tensors.push_back(t);
      }
    }
  };
}

#endif // DUNEURO_GMSH_TENSOR_READER_HH
