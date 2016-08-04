#ifndef DUNEURO_FUNCTION_HH
#define DUNEURO_FUNCTION_HH

#include <memory>

#include <dune/common/std/memory.hh>

namespace duneuro
{
  class Function
  {
    struct Base {
      virtual ~Base()
      {
      }
    };

    template <class T>
    struct Data : public Base {
      explicit Data(std::unique_ptr<T> d) : data(std::move(d))
      {
      }
      std::unique_ptr<T> data;
    };

  public:
    template <class T>
    explicit Function(std::unique_ptr<T> data)
        : data_(Dune::Std::make_unique<Data<T>>(std::move(data)))
    {
    }

    template <class T>
    T& cast()
    {
      Data<T>* ptr = dynamic_cast<Data<T>*>(data_.get());
      return *(ptr->data.get());
    }

    template <class T>
    const T& cast() const
    {
      const Data<T>* ptr = dynamic_cast<const Data<T>*>(data_.get());
      return *(ptr->data.get());
    }

  private:
    std::unique_ptr<Base> data_;
  };
}

#endif // DUNEURO_FUNCTION_HH
