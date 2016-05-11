#ifndef DUNEURO_IMAGE_HH
#define DUNEURO_IMAGE_HH

#include <memory>
#include <vector>

#include <duneuro/common/simple_structured_grid.hh>

namespace duneuro
{
  template <class T, int dim>
  class Image
  {
  public:
    using Grid = SimpleStructuredGrid<dim>;

    Image(std::shared_ptr<std::vector<T>> data, const Grid& grid) : data_(data), grid_(grid)
    {
      assert(data);
      if (data->size() != grid_.elements()) {
        DUNE_THROW(Dune::Exception, "image data size has to match the number of grid elements");
      }
    }

    const Grid& grid() const
    {
      return grid_;
    }

    const T& operator[](typename Grid::LinearIndex i) const
    {
      assert(i < data_->size());
      return (*data_)[i];
    }

    const T& operator[](const typename Grid::MultiIndex& i) const
    {
      auto li = grid_.toLinearIndex(i);
      assert(li < data_->size());
      return (*data_)[li];
    }

  private:
    std::shared_ptr<std::vector<T>> data_;
    Grid grid_;
  };
}

#endif // DUNEURO_IMAGE_HH
