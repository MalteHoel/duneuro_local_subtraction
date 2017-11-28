#ifndef DUNEURO_DIPOLE_STATISTICS_HH
#define DUNEURO_DIPOLE_STATISTICS_HH

#include <duneuro/common/dipole.hh>
#include <duneuro/common/fitted_driver_data.hh>
#include <duneuro/common/flags.hh>
#include <duneuro/common/volume_conductor_storage.hh>

namespace duneuro
{
  template <int dim>
  struct DipoleStatisticsInterface {
    using Coordinate = Dune::FieldVector<double, dim>;
    using DipoleType = Dipole<double, dim>;
    using TensorType = Dune::FieldMatrix<double, dim, dim>;

    virtual ~DipoleStatisticsInterface()
    {
    }

    virtual TensorType conductivity(const DipoleType& x) const = 0;
  };

  template <int dim, ElementType elementType, bool geometryAdaption = false>
  class FittedDipoleStatistics : public DipoleStatisticsInterface<dim>
  {
  public:
    using BaseT = DipoleStatisticsInterface<dim>;
    using Coordinate = typename BaseT::Coordinate;
    using DipoleType = typename BaseT::DipoleType;
    using TensorType = typename BaseT::TensorType;
    using VCStorage = VolumeConductorStorage<dim, elementType, geometryAdaption>;
    using VC = typename VCStorage::Type;
    using ElementSearch = KDTreeElementSearch<typename VC::GridView>;

    FittedDipoleStatistics(const FittedDriverData<dim>& data, const Dune::ParameterTree& config,
                           DataTree dataTree = DataTree())
        : config_(config)
        , volumeConductorStorage_(data, config.sub("volume_conductor"),
                                  dataTree.sub("volume_conductor"))
        , elementSearch_(std::make_shared<ElementSearch>(volumeConductorStorage_.get()->gridView()))
    {
    }

    virtual TensorType conductivity(const DipoleType& x) const override
    {
      return volumeConductorStorage_.get()->tensor(elementSearch_->findEntity(x.position()));
    }

  private:
    Dune::ParameterTree config_;
    VCStorage volumeConductorStorage_;
    std::shared_ptr<ElementSearch> elementSearch_;
  };
}

#endif // DUNEURO_DIPOLE_STATISTICS_HH
