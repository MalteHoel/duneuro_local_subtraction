#ifndef DUNEURO_FEATURE_MANAGER_HH
#define DUNEURO_FEATURE_MANAGER_HH

#include <sstream>

namespace duneuro {
struct FeatureCharacteristics {
  std::string description;
  std::vector<std::string> feature_papers;
};

class FeatureManager {
public:
  FeatureManager(const bool enable_experimental,
                 const Dune::ParameterTree &config)
      : enable_experimental_(enable_experimental),
        features_{
            {"duneuro",
             {"Please cite the following work related to duneuro:",
              {"Nüßing2018"}}},
            {"DUNE",
             {"Duneuro is based on and uses the numerical methods provided "
              "by DUNE. Here is a list of relevant publications:",
              {"Bastian2008", "Bastian2008b", "Bastian2010", "Blatt2007"}}},
            {"cg",
             {"Related to the Continuous Galerkin Finite Element Method "
              "(CG-FEM) to solve the EEG/MEG forward problems:",
              {"Vorwerk2014"}}},
            {"dg",
             {"Related to the Disontinuous Galerkin Finite Element Method "
              "(DG-FEM) to solve the EEG/MEG forward problems:",
              {"Engwer2017", "Piastra2018"}}},
            {"udg",
             {"Related to the Unfitted Disontinuous Galerkin Finite Element "
              "Method (UDG-FEM) to solve the EEG forward problems:",
              {"Nüßing2014"}}},
            {"cutfem", {"", {""}}},
            {"partial_integration",
             {"Related to the partial integration source model:",
              {"Bauer2014"}}},
            {"venant",
             {"Related to the St. Venant source model:", {"Bauer2014"}}},
            {"patch_based_venant", {"", {""}}},
            {"truncated_spatial_venant", {{""}}},
            {"subtraction",
             {"Related to the subtraction source model:",
              {"Wolters2007", "Drechsler2009"}}},
            {"localized_subtraction", {"", {""}}},
            {"whitney",
             {"Related to the Whitney source model:", {"Miinalainen2018"}}},
            {"transfer_matrix",
             {"Related to the transfer matrix approach to solve the EEG/MEG "
              "forward more efficiently:",
              {"Wolters2004"}}}},
        papers_{
            {"Bastian2008",
             "Bastian, P., Blatt, M., Dedner, A., Engwer, C., Klöfkorn, R., "
             "Ohlberger, M., Sander, O. (2008).\nA Generic Grid Interface for "
             "Parallel and Adaptive Scientific Computing. Part I: Abstract "
             "Framework.\nComputing, 82(2-3), pp. 103-119. (related to module: "
             "dune-grid)"},
            {"Bastian2008b",
             "Bastian, P., Blatt, M., Dedner, A., Engwer, C., Klöfkorn, R., "
             "Kornhuber, R., Ohlberger, M., Sander, O. (2008).\nA Generic "
             "Grid Interface for Parallel and Adaptive Scientific Computing. "
             "Part II: Implementation and Tests in DUNE.\nComputing, 82(2-3), "
             "pp. 121-138. (related to module: dune-grid)"},
            {"Bastian2010",
             "Bastian, P., Heimann, F., and Marnach, S. (2010).\nGeneric "
             "implementation of finite element methods in the Distributed and "
             "Unified Numerics Environment (DUNE).\nKybernetika, 46(2):294 "
             "315, (related to module: dune-pdelab)"},
            {"Blatt2007",
             "Blatt, M., Bastian, P. (2008).\nThe Iterative Solver Template "
             "Library.\nIn B. Kåström, E. Elmroth, J. Dongarra and J. "
             "Wasniewski, Applied Parallel Computing. State of the Art in "
             "Scientific Computing.\nVolume 4699 of Lecture Notes in "
             "Scientific Computing, pp. 666-675. Springer, (related to module: "
             "dune-istl)"},
            {"Bauer2014",
             "Bauer, M., Pursiainen, S., Vorwerk, J., Köstler, H., Wolters, "
             "C.H. (2015).\nComparison Study for Whitney (Raviart-Thomas) "
             "Type Source Models in Finite Element Method Based EEG Forward "
             "Modeling.\nIEEE Trans. Biomed. Eng.,62(11):2648–56"},
            {"Drechsler2009",
             "Drechsler, F., Wolters, C. H., Dierkes, T., Si, H., and "
             "Grasedyck, L. (2009).\nA full subtraction approach for finite "
             "element method based source analysis using constrained Delaunay "
             "tetrahedralisation.\nNeuroImage 46, 1055–1065."},
            {"Engwer2017",
             "Engwer, C., Vorwerk, J., Ludewig, J., and Wolters, C. H. "
             "(2017).\nA discontinuous galerkin method to solve the EEG "
             "forward problem using the subtraction approach.\nSIAM J. Sci. "
             "Comput. 39, B138–B164. (EEG)"},
            {"Miinalainen2018",
             "Miinalainen, T., Rezaei, A., Defne, U., Nüssing, A. & Engwer, "
             "C. & Wolters, C. H., Pursiainen, S. (2019).\nA realistic, "
             "accurate and fast source modeling approach for the EEG forward "
             "problem.\nNeuroimage 184:56-67."},
            {"Nüßing2016",
             "Nüßing, A., Wolters, C. H., Brinck, H., and Engwer, C. "
             "(2016).\nThe unfitted discontinuous Galerkin method for solving "
             "the EEG forward problem.\nIEEE Trans. Biomed. Eng. 63, "
             "2564–2575."},
            {"Nüßing2018",
             "Nüßing, A., Piastra, M.C., Schrader, S., Miinalainen, T., "
             "Pursiainen, S., Brinck, H., Wolters, C.H., Engwer, C. "
             "(2019).\nduneuro - A software toolbox for forward modeling in "
             "neuroscience.\narXiv, eprint={1901.02874}"},
            {"Piastra2018",
             "Piastra, M. C., Nüßing, A., Vorwerk, J., Bornfleth, H., "
             "Oostenveld, R., Engwer, C., Wolters, C. H. (2018).\nThe "
             "Discontinuous Galerkin Finite Element Method for Solving the "
             "MEG and the Combined MEG/EEG Forward Problem.\nFrontiers in "
             "Neuroscience 12, 30"},
            {"Vorwerk2014",
             "Vorwerk, J., Cho, J.-H., Rampp, S., Hamer, H., Knösche,T.R., "
             "Wolters,C.H. (2014).\nA Guideline for Head Volume Conductor "
             "Modeling in EEG and MEG.\nNeuroImage, 100, pp.590-607. (MEG)"},
            {"Wolters2004",
             "Wolters, C. H., Grasedyck, L., and Hackbusch, W. "
             "(2004).\nEfficient computation of lead field bases and "
             "influence matrix for the FEM-based EEG and MEG inverse "
             "problem.\nInverse Prob. 20:1099"},
            {"Wolters2007",
             "Wolters, C.H., Köstler, H., Möller, C., Härdtlein, J., "
             "Grasedyck, L., Hackbusch, W.(2007).\nNumerical mathematics of "
             "the subtraction approach for the modeling of a current dipole "
             "in EEG source reconstruction using finite element head "
             "models.\nSIAM J. on Scientific Computing, 30(1):24-45."}},
        experimental_features_{"cutfem", "patch_based_venant",
                               "truncated_spatial_venant",
                               "localized_subtraction"} {
    check_feature(config);
  }

  void check_feature(const Dune::ParameterTree &config) {
    std::string feature_name;
    if (config.hasKey("solver_type")) {
      feature_name = config.get<std::string>("solver_type");
    } else if (config.hasSub("source_model")) {
      feature_name = config.get<std::string>("source_model.type");
    }
    check_feature(feature_name);
  }

  void check_feature(const std::string &feature_name) {
    if (!feature_name.empty()) {
      check_access(feature_name);
      relevant_features_.push_back(feature_name);
    }
  }

  void check_access(const std::string feature_name) const {
    if (get_experimental(feature_name)) {
      if (!enable_experimental_) {
        DUNE_THROW(Dune::Exception, "The requested feature "
                                        << feature_name << "is experimental.");
      } else {
        std::cout << "WARNING: The requested feature " << feature_name
                  << "is experimental." << std::endl;
      }
    }
  }

  void print_citations() {
    relevant_features_.sort();
    relevant_features_.unique();
    relevant_features_.insert(relevant_features_.begin(), {"duneuro", "DUNE"});
    for (const std::string &feature_id : relevant_features_) {
      std::cout << get_description(feature_id);
      std::cout << "\n\n";
      std::vector<std::string> feature_papers = get_feature_papers(feature_id);
      for (const std::string &paper_id : feature_papers) {
        std::cout << get_citation(paper_id);
        std::cout << "\n\n";
      }
    }
  }

  bool get_experimental(const std::string &feature_name) const {
    return (experimental_features_.find(feature_name) !=
            experimental_features_.end());
  }

  std::vector<std::string>
  get_feature_papers(const std::string &feature_name) const {
    return features_.at(feature_name).feature_papers;
  }

  std::string get_description(const std::string &feature_name) const {
    return features_.at(feature_name).description;
  }

  std::string get_citation(const std::string &paper_id) const {
    return papers_.at(paper_id);
  }

private:
  const bool enable_experimental_;
  const std::map<std::string, FeatureCharacteristics> features_;
  const std::map<std::string, std::string> papers_;
  std::list<std::string> relevant_features_;
  std::set<std::string> experimental_features_;
};
} // namespace duneuro

#endif // DUNEURO_FEATURE_MANAGER_HH
