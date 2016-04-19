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
    }

    const Grid& grid() const
    {
      return grid_;
    }

    const T& operator[](typename Grid::LinearIndex i) const
    {
      return (*data_)[i];
    }

    const T& operator[](const typename Grid::MultiIndex& i) const
    {
      return (*data_)[grid_.toLinearIndex(i)];
    }

  private:
    std::shared_ptr<std::vector<T>> data_;
    Grid grid_;
  };
}

#endif // DUNEURO_IMAGE_HH
