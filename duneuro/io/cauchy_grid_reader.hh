// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_CAUCHYGRIDREADER_HH
#define DUNEURO_CAUCHYGRIDREADER_HH

#include <fstream>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/gridfactory.hh>

#include <duneuro/io/cauchy_utilities.hh>

namespace duneuro
{
  namespace CauchyDetail
  {
    template <class Grid>
    class CauchyGridParser
    {
    public:
      CauchyGridParser(Dune::GridFactory<Grid>& factory, bool verbose)
          : factory_(factory), verbose_(verbose), knotenAnzahl_(0), elementAnzahl_(0)
      {
      }
      void read(const std::string& filename);

    private:
      void readGeometryFile(std::istream& stream);
      void readSteuerkarte(std::istream& stream);
      void readKoordinatenkarte(std::istream& stream);
      void readElementKnotenkarte(std::istream& stream);
      Dune::GridFactory<Grid>& factory_;
      bool verbose_;
      unsigned int knotenAnzahl_;
      unsigned int elementAnzahl_;
    };

    template <class Grid>
    void CauchyGridParser<Grid>::read(const std::string& filename)
    {
      std::ifstream stream(filename);
      if (!stream) {
        DUNE_THROW(Dune::IOError, "file " << filename << " could not be opended");
      }
      readGeometryFile(stream);
    }

    template <class Grid>
    void CauchyGridParser<Grid>::readGeometryFile(std::istream& stream)
    {
      readLine(stream, "BOI - GEOMETRIEFILE");
      ignoreSeparatorLine(stream, 2);
      readSteuerkarte(stream);
      ignoreSeparatorLine(stream, 2);
      readKoordinatenkarte(stream);
      ignoreSeparatorLine(stream, 2);
      readElementKnotenkarte(stream);
      ignoreSeparatorLine(stream, 2);
      readLine(stream, "EOI - GEOMETRIEFILE");
    }

    template <class Grid>
    void CauchyGridParser<Grid>::readSteuerkarte(std::istream& stream)
    {
      readLine(stream, "BOI - STEUERKARTE");
      for (;;) {
        std::string line;
        std::getline(stream, line, '\n');
        if (line == "EOI - STEUERKARTE") {
          break;
        }
        // get key
        std::stringstream sstr;
        sstr << line;
        std::getline(sstr, line, ':');
        // strip line
        strip(line);
        if (line == "ANZAHL DER KNOTEN") {
          sstr >> knotenAnzahl_;
        } else if (line == "ANZAHL DER ELEMENTE") {
          sstr >> elementAnzahl_;
        }
      }
    }

    template <class Grid>
    void CauchyGridParser<Grid>::readKoordinatenkarte(std::istream& stream)
    {
      readLine(stream, "BOI - KOORDINATENKARTE");
      for (unsigned int i = 0; i < knotenAnzahl_; ++i) {
        Dune::FieldVector<double, Grid::dimension> v;
        for (unsigned int j = 0; j < Grid::dimension; ++j) {
          stream >> v[j];
        }
        factory_.insertVertex(v);
      }
      std::string line;
      std::getline(stream, line, '\n');
      // stream.ignore(std::numeric_limits<std::size_t>::max(),'\n');
      readLine(stream, "EOI - KOORDINATENKARTE");
    }

    template <class Grid>
    void CauchyGridParser<Grid>::readElementKnotenkarte(std::istream& stream)
    {
      readLine(stream, "BOI - ELEMENTKNOTENKARTE");
      for (unsigned int i = 0; i < elementAnzahl_; ++i) {
        // retrieve element type
        std::string str;
        std::getline(stream, str, ':');
        std::stringstream convStream(str);
        int type;
        convStream >> type;
        Dune::GeometryType gtype;
        switch (type) {
        case 303: gtype = Dune::GeometryTypes::simplex(Grid::dimension); break;
        case 323: gtype = Dune::GeometryTypes::cube(Grid::dimension); break;
        default: DUNE_THROW(Dune::IOError, "unknown geometry type: " << type);
        }
        // get element indices
        std::getline(stream, str, '\n');
        // remove first space
        str.erase(0, 1);
        std::vector<unsigned int> indices;
        unsigned int v;
        unsigned int index = 0;
        while (index < str.size()) {
          // retrieve substring, 6 characters for FORTRAN
          std::string currentString = str.substr(index, 6);
          strip(currentString);
          std::stringstream currentIndex(currentString);
          currentIndex >> v;
          // v-1 because FORTRAN indices start with 1!
          indices.push_back(v - 1);
          index += 6;
        }
        factory_.insertElement(gtype, indices);
      }
      readLine(stream, "EOI - ELEMENTKNOTENKARTE");
    }
  }
  /** \brief creates grid from a gives SimBio-Cauchy field format */
  template <class Grid>
  class CauchyGridReader
  {
  public:
    static void read(Dune::GridFactory<Grid>& factory, const std::string& filename,
                     bool verbose = true);
  };

  template <class Grid>
  void CauchyGridReader<Grid>::read(Dune::GridFactory<Grid>& factory, const std::string& filename,
                                    bool verbose)
  {
    CauchyDetail::CauchyGridParser<Grid> parser(factory, verbose);
    parser.read(filename);
  }
}
#endif // DUNEURO_CAUCHYGRIDREADER_HH
