// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_NIFTI_IMAGE_WRITER_HH
#define DUNEURO_NIFTI_IMAGE_WRITER_HH

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <duneuro/io/nifti_utilities.hh>

namespace duneuro
{
  class NiftiImageWriter
  {
  public:
    template <int dim>
    static void write(const std::string& filename, const std::array<unsigned int, dim>& cells,
                      const std::vector<char>& data)
    {
      int dims[8];
      int datatype = DT_INT8;
      std::fill(dims, dims + 8, 0);
      dims[0] = dim;
      std::copy(cells.begin(), cells.end(), dims + 1);
      std::shared_ptr<nifti_image> image(nifti_make_new_nim(dims, datatype, 1));
      // use only .nii file
      image->nifti_type = 1;
      std::copy(data.begin(), data.end(), reinterpret_cast<char*>(image->data));
      nifti_set_filenames(image.get(), filename.c_str(), 1, 1);
      nifti_image_write(image.get());
    }
  };
}

#endif // DUNEURO_NIFTI_IMAGE_WRITER_HH
