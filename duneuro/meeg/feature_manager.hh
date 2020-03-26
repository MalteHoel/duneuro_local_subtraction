#ifndef DUNEURO_FEATURE_MANAGER_HH
#define DUNEURO_FEATURE_MANAGER_HH

#include <sstream>

namespace duneuro {
struct FeatureCharacteristics {
  bool experimental;
  std::list<std::string> feature_papers;
};

struct PaperCharacteristics {
  std::string citation;
  std::string comment;
};

class FeatureManager {
public:
  FeatureManager(const bool enable_experimental,
                 const Dune::ParameterTree &config)
      : enable_experimental_(enable_experimental),
        features_{{"cg", {false, {"Vorwerk2014"}}},
                  {"dg", {false, {"Engwer2017", "Piastra2018"}}},
                  {"udg", {false, {"Nüßing2014"}}},
                  {"cutfem", {true, {""}}},
                  {"partial_integration", {false, {"Bauer2014"}}},
                  {"venant", {false, {"Bauer2014"}}},
                  {"patch_based_venant", {true, {""}}},
                  {"truncated_spatial_venant", {true, {""}}},
                  {"subtraction", {false, {"Wolters2007", "Drechsler2009"}}},
                  {"localized_subtraction", {true, {""}}},
                  {"whitney", {false, {"Miinalainen2018"}}},
                  {"transfer_matrix", {false, {"Wolters2004"}}}},
        papers_{
            {"Bastian2008",
             {"Bastian, P., Blatt, M., Dedner, A., Engwer, C., Klöfkorn, R., "
              "Ohlberger, M., Sander, O. (2008).\nA Generic Grid Interface for "
              "Parallel and Adaptive Scientific Computing. Part I: Abstract "
              "Framework.\nComputing, 82(2-3), pp. 103-119.",
              "dune-grid"}},
            {"Bastian2008b",
             {"Bastian, P., Blatt, M., Dedner, A., Engwer, C., Klöfkorn, R., "
              "Kornhuber, R., Ohlberger, M., Sander, O. (2008).\nA Generic "
              "Grid Interface for Parallel and Adaptive Scientific Computing. "
              "Part II: Implementation and Tests in DUNE.\nComputing, 82(2-3), "
              "pp. 121-138.",
              "dune-grid"}},
            {"Bastian2010",
             {"Bastian, P., Heimann, F., and Marnach, S. (2010).\nGeneric "
              "implementation of finite element methods in the Distributed and "
              "Unified Numerics Environment (DUNE).\nKybernetika, 46(2):294 "
              "315",
              "dune-pdelab"}},
            {"Blatt2007",
             {"Blatt, M., Bastian, P. (2008).\nThe Iterative Solver Template "
              "Library.\nIn B. Kåström, E. Elmroth, J. Dongarra and J. "
              "Wasniewski, Applied Parallel Computing. State of the Art in "
              "Scientific Computing.\nVolume 4699 of Lecture Notes in "
              "Scientific Computing, pp. 666-675. Springer",
              "dune-istl"}},
            {"Bauer2014",
             {"Bauer, M., Pursiainen, S., Vorwerk, J., Köstler, H., Wolters, "
              "C.H. (2015).\nComparison Study for Whitney (Raviart-Thomas) "
              "Type Source Models in Finite Element Method Based EEG Forward "
              "Modeling.\nIEEE Trans. Biomed. Eng.,62(11):2648–56",
              ""}},
            {"Drechsler2009",
             {"Drechsler, F., Wolters, C. H., Dierkes, T., Si, H., and "
              "Grasedyck, L. (2009).\nA full subtraction approach for finite "
              "element method based source analysis using constrained Delaunay "
              "tetrahedralisation.\nNeuroImage 46, 1055–1065.",
              ""}},
            {"Engwer2017",
             {"Engwer, C., Vorwerk, J., Ludewig, J., and Wolters, C. H. "
              "(2017).\nA discontinuous galerkin method to solve the EEG "
              "forward problem using the subtraction approach.\nSIAM J. Sci. "
              "Comput. 39, B138–B164.",
              "EEG"}},
            {"Miinalainen2018",
             {"Miinalainen, T., Rezaei, A., Defne, U., Nüssing, A. & Engwer, "
              "C. & Wolters, C. H., Pursiainen, S. (2019).\nA realistic, "
              "accurate and fast source modeling approach for the EEG forward "
              "problem.\nNeuroimage 184:56-67.",
              ""}},
            {"Nüßing2016",
             {"Nüßing, A., Wolters, C. H., Brinck, H., and Engwer, C. "
              "(2016).\nThe unfitted discontinuous Galerkin method for solving "
              "the EEG forward problem.\nIEEE Trans. Biomed. Eng. 63, "
              "2564–2575.",
              ""}},
            {"Nüßing2018",
             {"Nüßing, A., Piastra, M.C., Schrader, S., Miinalainen, T., "
              "Pursiainen, S., Brinck, H., Wolters, C.H., Engwer, C. "
              "(2019).\nduneuro - A software toolbox for forward modeling in "
              "neuroscience.\narXiv, eprint={1901.02874}",
              ""}},
            {"Piastra2018",
             {"Piastra, M. C., Nüßing, A., Vorwerk, J., Bornfleth, H., "
              "Oostenveld, R., Engwer, C., Wolters, C. H. (2018).\nThe "
              "Discontinuous Galerkin Finite Element Method for Solving the "
              "MEG and the Combined MEG/EEG Forward Problem.\nFrontiers in "
              "Neuroscience 12, 30",
              ""}},
            {"Vorwerk2014",
             {"Vorwerk, J., Cho, J.-H., Rampp, S., Hamer, H., Knösche,T.R., "
              "Wolters,C.H. (2014).\nA Guideline for Head Volume Conductor "
              "Modeling in EEG and MEG.\nNeuroImage, 100, pp.590-607",
              "MEG"}},
            {"Wolters2004",
             {"Wolters, C. H., Grasedyck, L., and Hackbusch, W. "
              "(2004).\nEfficient computation of lead field bases and "
              "influence matrix for the FEM-based EEG and MEG inverse "
              "problem.\nInverse Prob. 20:1099",
              ""}},
            {"Wolters2007",
             {"Wolters, C.H., Köstler, H., Möller, C., Härdtlein, J., "
              "Grasedyck, L., Hackbusch, W.(2007).\nNumerical mathematics of "
              "the subtraction approach for the modeling of a current dipole "
              "in EEG source reconstruction using finite element head "
              "models.\nSIAM J. on Scientific Computing, 30(1):24-45.",
              ""}}} {
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
      add_citation(feature_name);
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

  void add_citation(const std::string feature_name) {
    for (const std::string &feature_paper : get_feature_papers(feature_name)) {
      relevant_feature_papers_.push_back(feature_paper);
    }
  }

  void print_citations() {
    std::cout << "Please cite the following work related to duneuro: "
              << std::endl;
    print_paper_information("Nüßing2018");
    std::cout << "Duneuro is a module of dune, here is a list of relevant "
                 "publications: "
              << std::endl;
    for (const std::string &paper_id :
         {"Bastian2008", "Bastian2008b", "Bastian2010", "Blatt2007"}) {
      print_paper_information(paper_id);
    }
    std::cout
        << "Here is a list of publications related to the features (finite "
           "element discretizations, source models,...) that you used: "
        << std::endl;
    relevant_feature_papers_.sort();
    relevant_feature_papers_.unique();
    for (const std::string &paper_id : relevant_feature_papers_) {
      print_paper_information(paper_id);
    }
  }

  void print_paper_information(const std::string &paper_id) const {
    std::cout << get_citation(paper_id);
    if (!get_comment(paper_id).empty()) {
      std::cout << " (" << get_comment(paper_id) << ") ";
    }
    std::cout << std::endl;
  }

  bool get_experimental(const std::string &feature_name) const {
    return features_.at(feature_name).experimental;
  }

  std::list<std::string> get_feature_papers(const std::string &feature_name) const {
    return features_.at(feature_name).feature_papers;
  }

  std::string get_citation(const std::string &paper_id) const {
    return papers_.at(paper_id).citation;
  }

  std::string get_comment(const std::string &paper_id) const {
    return papers_.at(paper_id).comment;
  }

private:
  const bool enable_experimental_;
  const std::map<std::string, FeatureCharacteristics> features_;
  const std::map<std::string, PaperCharacteristics> papers_;
  std::list<std::string> relevant_feature_papers_;
};
} // namespace duneuro

#endif // DUNEURO_FEATURE_MANAGER_HH
