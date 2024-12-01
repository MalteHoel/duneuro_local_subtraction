// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_GMSH_WRITER_HH
#define DUNEURO_GMSH_WRITER_HH

#include <memory>

#include <dune/grid/io/file/gmshwriter.hh>

namespace duneuro
{
  class GmshWriter
  {
  public:
    template <class G>
    static void write(const std::string& filename, const G& grid, const std::vector<char>& labels)
    {
      typename G::LeafGridView gv(grid.leafGridView());
      Dune::GmshWriter<typename G::LeafGridView> writer(gv);
      std::vector<int> il(labels.begin(), labels.end());
      writer.write(filename, il);
    }

  private:
  };
}

#endif DUNEURO_GMSH_WRITER_HH
