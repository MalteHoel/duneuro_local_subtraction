// SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
#ifndef DUNEURO_NIFTI_UTILITIES_HH
#define DUNEURO_NIFTI_UTILITIES_HH

#include <nifti1_io.h>

namespace duneuro
{
  struct NiftiDeleter {
    void operator()(nifti_image* img)
    {
      if (img) {
        nifti_image_free(img);
      }
    }
  };
}

#endif // DUNEURO_NIFTI_UTILITIES_HH
