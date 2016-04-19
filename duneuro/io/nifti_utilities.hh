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
