#ifndef DUNEURO_FUNCTION_HH
#define DUNEURO_FUNCTION_HH

#include <memory>

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
      explicit Data(std::shared_ptr<T> d) : data(d)
      {
      }
      std::shared_ptr<T> data;
    };

  public:
    template <class T>
    explicit Function(std::shared_ptr<T> data) : data_(new Data<T>(data))
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
    std::shared_ptr<Base> data_;
  };
}

#endif // DUNEURO_FUNCTION_HH
