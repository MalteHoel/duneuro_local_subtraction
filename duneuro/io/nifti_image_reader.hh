#ifndef DUNEURO_NIFTI_IMAGE_READER_HH
#define DUNEURO_NIFTI_IMAGE_READER_HH

#include <nifti1_io.h>

#include <dune/common/exceptions.hh>

#include <duneuro/common/image.hh>
#include <duneuro/common/simple_structured_grid.hh>
#include <duneuro/io/nifti_utilities.hh>

namespace duneuro
{
  class NiftiImageReader
  {
  public:
    template <int dim, class I>
    static void read(const std::string& filename, I out, std::array<unsigned int, dim>& cells)
    {
      std::shared_ptr<nifti_image> nim(nifti_image_read(filename.c_str(), true), NiftiDeleter());
      // nifti_image_infodump(nim.get());
      std::copy(nim->dim + 1, nim->dim + dim + 1, cells.begin());
      if (nim->datatype == NIFTI_TYPE_FLOAT32) {
        std::copy(static_cast<float*>(nim->data), static_cast<float*>(nim->data) + nim->nvox, out);
      } else if (nim->datatype == NIFTI_TYPE_FLOAT64) {
        std::copy(static_cast<double*>(nim->data), static_cast<double*>(nim->data) + nim->nvox,
                  out);
      } else if (nim->datatype == NIFTI_TYPE_INT8) {
        std::copy(static_cast<signed char*>(nim->data),
                  static_cast<signed char*>(nim->data) + nim->nvox, out);
      } else {
        DUNE_THROW(Dune::Exception, "unknown nifti data type");
      }
    }

    template <class T, int dim>
    static std::shared_ptr<Image<T, dim>> read(const std::string& filename)
    {
      auto data = std::make_shared<std::vector<T>>();
      std::array<unsigned int, dim> cells;
      read<dim>(filename, std::back_inserter(*data), cells);
      return std::make_shared<Image<T, dim>>(data, SimpleStructuredGrid<dim>(cells));
    }
  };
}

#endif // DUNEURO_NIFTI_IMAGE_READER_HH
