#ifndef DUNEURO_FEATURE_MANAGER_HH
#define DUNEURO_FEATURE_MANAGER_HH

#include <sstream>

namespace duneuro {
struct FeatureCharacteristics {
  bool experimental;
  std::string citation;
};

class FeatureManager {
public:
  FeatureManager(const bool enable_experimental,
                 const Dune::ParameterTree &config)
      : enable_experimental_(enable_experimental),
        features_{{"cg", {false, ""}},
                  {"dg", {false, "X, Y"}},
                  {"udg", {false, "X"}},
                  {"cutfem", {true, " "}},
                  {"partial_integration", {false, ""}},
                  {"venant", {false, ""}},
                  {"pathc_based_venant", {true, ""}},
                  {"truncated_spatial_venant", {true, ""}},
                  {"subtraction", {false, ""}},
                  {"localized_subtraction", {true, ""}},
                  {"whitney", {false, "X"}}} {
    check_feature(config);
  }

  void check_feature(const Dune::ParameterTree &config) {
    std::string key;
    if (config.hasKey("solver_type")) {
      key = config.get<std::string>("solver_type");
    } else if (config.hasSub("source_model")) {
      key = config.get<std::string>("source_model.type");
    }
    if (!key.empty()) {
      check_access(key);
      add_citation(key);
    }
  }

  void check_access(const std::string key) const {
    if (get_experimental(key)) {
      if (!enable_experimental_) {
        DUNE_THROW(Dune::Exception,
                   "The requested feature " << key << "is experimental.");
      } else {
        std::cout << "WARNING: The requested feature " << key
                  << "is experimental." << std::endl;
      }
    }
  }

  void add_citation(const std::string key) {
    features_used.push_back(get_citation(key));
  }

  void print_citations() const {
    std::cout << "Please cite the following works related to the features "
                 "you used: "
              << std::endl;
    for (const std::string &citation : features_used) {
      std::cout << citation << std::endl;
    }
  }

  bool get_experimental(std::string key) const {
    return features_.at(key).experimental;
  }

  std::string get_citation(std::string key) const {
    return features_.at(key).citation;
  }

private:
  const std::map<std::string, FeatureCharacteristics> features_;
  const bool enable_experimental_;
  std::list<std::string> features_used;
};
} // namespace duneuro

#endif // DUNEURO_FEATURE_MANAGER_HH
