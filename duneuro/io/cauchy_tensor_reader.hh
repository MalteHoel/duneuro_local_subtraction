#ifndef DUNEURO_CAUCHYTENSORREADER_HH
#define DUNEURO_CAUCHYTENSORREADER_HH

#include <fstream>
#include <iostream>

#include <dune/common/fmatrix.hh>

#include <duneuro/io/cauchy_utilities.hh>

namespace duneuro
{
  namespace CauchyDetail
  {
    template <class T>
    struct Cauchy3DTensor {
      typedef Dune::FieldMatrix<double, 3, 3> MatrixType;
      MatrixType createMatrix() const
      {
        MatrixType matrix;
        // store six tensor values (xx,yy,zz,xy,yz,xz)
        matrix[0][0] = values[0]/1000;
        matrix[1][1] = values[1]/1000;
        matrix[2][2] = values[2/1000];
        matrix[0][1] = matrix[1][0] = values[3]/1000;
        matrix[1][2] = matrix[2][1] = values[4]/1000;
        matrix[0][2] = matrix[2][0] = values[5]/1000;
        return matrix;
      }
      std::array<T, 6> values;
    };

    template <class T>
    struct Cauchy2DTensor {
      typedef Dune::FieldMatrix<double, 2, 2> MatrixType;
      MatrixType createMatrix() const
      {
        MatrixType matrix;
        // store six tensor values (xx,yy,xy)
        matrix[0][0] = values[0]/1000;
        matrix[1][1] = values[1]/1000;
        matrix[0][1] = matrix[1][0] = values[2]/1000;
        return matrix;
      }
      std::array<T, 3> values;
    };

    template <class T, int dim>
    struct SelectCauchyTensor;

    template <class T>
    struct SelectCauchyTensor<T, 2> {
      using Type = Cauchy2DTensor<T>;
    };

    template <class T>
    struct SelectCauchyTensor<T, 3> {
      using Type = Cauchy3DTensor<T>;
    };

    template <class T>
    std::istream& operator>>(std::istream& stream, Cauchy3DTensor<T>& tensor)
    {
      for (std::size_t i = 0; i < tensor.values.size(); ++i) {
        stream >> tensor.values[i];
      }
      return stream;
    }

    template <class T>
    std::istream& operator>>(std::istream& stream, Cauchy2DTensor<T>& tensor)
    {
      for (std::size_t i = 0; i < tensor.values.size(); ++i) {
        stream >> tensor.values[i];
      }
      return stream;
    }

    template <class Grid>
    class CauchyTensorParser
    {
    public:
      enum { dim = Grid::dimension };
      typedef typename Grid::ctype ctype;
      typedef typename SelectCauchyTensor<ctype, dim>::Type CauchyTensor;
      typedef typename CauchyTensor::MatrixType TensorMatrixType;

      explicit CauchyTensorParser(bool verbose) : verbose_(verbose)
      {
      }
      template <class I>
      void read(const std::string& filename, I output);

    private:
      template <class I>
      void readTensorValueFile(std::istream& stream, I output);
      template <class I>
      void readTensor(std::istream& stream, I output);
      bool verbose_;
    };

    template <class Grid>
    template <class I>
    void CauchyTensorParser<Grid>::read(const std::string& filename, I output)
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "file " << filename << " could not be opended");
      }
      readTensorValueFile(stream, output);
    }

    template <class Grid>
    template <class I>
    void CauchyTensorParser<Grid>::readTensorValueFile(std::istream& stream, I output)
    {
      readLine(stream, "BOI - TENSORVALUEFILE");
      ignoreSeparatorLine(stream, 2);
      readTensor(stream, output);
      ignoreSeparatorLine(stream, 2);
      readLine(stream, "EOI - TENSORVALUEFILE");
    }

    template <class Grid>
    template <class I>
    void CauchyTensorParser<Grid>::readTensor(std::istream& stream, I output)
    {
      readLine(stream, "BOI - TENSOR");
      for (unsigned int count = 1;; ++count) {
        std::string line;
        std::getline(stream, line, '\n');
        if (line == "EOI - TENSOR") {
          break;
        }
        // now in a tensor line
        // extract second line
        std::string secondLine;
        std::getline(stream, secondLine, '\n');
        std::stringstream sstr;
        sstr << line << secondLine;
        // retrieve index
        unsigned int index;
        sstr >> index;
        if (index != count) {
          DUNE_THROW(Dune::IOError, "indices in tensor file not consecutive (found "
                                        << index << " but expected " << count << ")");
        }
        CauchyTensor tensor;
        sstr >> tensor;
        *output++ = tensor.createMatrix();
      }
    }
  }

  template <class Grid>
  class CauchyTensorReader
  {
  public:
    enum { dimension = Grid::dimension };
    typedef CauchyDetail::CauchyTensorParser<Grid> Parser;

    template <class I>
    static void read(const std::string& filename, I output, bool verbose = true);
  };

  template <class Grid>
  template <class I>
  void CauchyTensorReader<Grid>::read(const std::string& filename, I output, bool verbose)
  {
    Parser parser(verbose);
    parser.read(filename, output);
  }
}

#endif // DUNEURO_CAUCHYTENSORREADER_HH
