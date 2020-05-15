#ifndef DUNEURO_FEATURE_MANAGER_HH
#define DUNEURO_FEATURE_MANAGER_HH

#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <string>

namespace duneuro {
struct FeatureCharacteristics {
  std::string description;
  std::vector<std::string> feature_papers;
};

class FeatureManager {
public:
  FeatureManager(const bool enable_experimental)
      : enable_experimental_(enable_experimental),
        features_{
            {"duneuro",
             {"Please cite the following work related to duneuro:",
              {"Nüßing2019"}}},
            {"DUNE",
             {"Duneuro is based on and uses the numerical methods provided "
              "by DUNE. Here is a list of relevant publications:",
              {"Bastian2008", "Bastian2008b", "Bastian2008c", "Blatt2007"}}},
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
            {"cutfem", {"", {}}},
            {"partial_integration",
             {"Related to the partial integration source model:",
              {"Bauer2015"}}},
            {"venant",
             {"Related to the St. Venant source model:", {"Bauer2015"}}},
            {"patch_based_venant", {"", {}}},
            {"truncated_spatial_venant", {"", {}}},
            {"subtraction",
             {"Related to the subtraction source model:",
              {"Wolters2007", "Drechsler2009"}}},
            {"localized_subtraction", {"", {}}},
            {"whitney",
             {"Related to the Whitney source model:", {"Miinalainen2018"}}},
            {"transfer_matrix",
             {"Related to the transfer matrix approach to solve the EEG/MEG "
              "forward more efficiently:",
              {"Wolters2004"}}}},
        papers_{
            {"Bastian2008",
             "P. Bastian, M. Blatt, A. Dedner, C. Engwer, R. Klöfkorn, M. "
             "Ohlberger, and O. Sander (2008).\nA Generic Grid Interface for "
             "Parallel and Adaptive Scientific Computing. Part I: Abstract "
             "Framework.\nComputing, 82(2-3), pp. 103-119."},
            {"Bastian2008b",
             "P. Bastian, M. Blatt, A. Dedner, C. Engwer, R. Klöfkorn, R. "
             "Kornhuber, M. Ohlberger, and O. Sander (2008).\nA Generic Grid "
             "Interface for Parallel and Adaptive Scientific Computing. Part "
             "II: Implementation and Tests in DUNE.\nComputing, 82(2-3), pp. "
             "121-138."},
            {"Bastian2008c",
             "P. Bastian, M. Blatt (2008). On the Generic Parallelisation of "
             "Iterative Solvers for the Finite Element Method In Int. J. "
             "Computational Science and Engineering,4(1), pp. 56-69."},
            {"Blatt2007",
             "M. Blatt, P. Bastian (2007).\nThe Iterative Solver Template "
             "Library.\nIn B. Kåström, E. Elmroth, J. Dongarra and J. "
             "Wasniewski, Applied Parallel Computing. State of the Art in "
             "Scientific Computing.\nVolume 4699 of Lecture Notes in "
             "Scientific Computing, pp. 666-675. Springer."},
            {"Bauer2015",
            "M. Bauer, S. Pursiainen, J. Vorwerk, H. Köstler, "
            "and C. H. Wolters (2015).\nComparison Study for "
            "Whitney (Raviart–Thomas)-Type Source Models in "
            "Finite-Element-Method-Based EEG Forward "
            "Modeling.\nIEEE Trans. Biomed. Eng., 62(11), pp. "
            "2648-2656, 2015, doi: 10.1109/TBME.2015.2439282."},
            {"Drechsler2009",
             "F. Drechsler, C.H. Wolters, T. Dierkes, H. Si, and L. Grasedyck "
             "(2009).\nA full subtraction approach for finite element method "
             "based source analysis using constrained Delaunay "
             "tetrahedralisation.\nNeuroImage 46(4), pp. 1055–1065, "
             "doi:10.1016/j.neuroimage.2009.02.024."},
            {"Engwer2017",
             "C. Engwer, J. Vorwerk, J. Ludewig, and C. H. Wolters (2017).\nA "
             "discontinuous Galerkin method to solve the EEG forward problem "
             "using the subtraction approach.\nSIAM J. Sci. Comput. 39(1), "
             "B138–B164, doi: 10.1137/15M1048392."},
            {"Miinalainen2018",
             "T. Miinalainen, A. Rezaei, D. Us, A. Nüßing, C. Engwer, C. H. "
             "Wolters, and S. Pursiainen (2019).\nA realistic, accurate and "
             "fast source modeling approach for the EEG forward "
             "problem.\nNeuroimage 184, pp. 56-67, doi: "
             "10.1016/j.neuroimage.2018.08.054."},
            {"Nüßing2016",
             "A. Nüßing, C. H. Wolters, H. Brinck, and C. Engwer (2016).\nThe "
             "unfitted discontinuous Galerkin method for solving the EEG "
             "forward problem.\nIEEE Trans. Biomed. Eng. 63(12), pp. "
             "2564–2575, doi:  10.1109/TBME.2016.2590740."},
            {"Nüßing2019",
             "A. Nüßing, M.C. Piastra, S. Schrader, T. Miinalainen, S. "
             "Pursiainen, H. Brinck, C. H. Wolters, and C. Engwer "
             "(2019).\nduneuro - A software toolbox for forward modeling in "
             "neuroscience.\narXiv, eprint={1901.02874}"},
            {"Piastra2018",
             "M. C. Piastra, A. Nüßing, J. Vorwerk, H. Bornfleth, R. "
             "Oostenveld, C. Engwer, and C. H. Wolters (2018).\nThe "
             "Discontinuous Galerkin Finite Element Method for Solving the MEG "
             "and the Combined MEG/EEG Forward Problem.\nFrontiers in "
             "Neuroscience 12, 30, doi: 10.3389/fnins.2018.00030."},
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
    if (enable_experimental) {
      std::cout << "WARNING: You have enabled to use untested features which "
                   "may lead to unreliable results."
                << std::endl;
    }
  }

  void check_feature(Dune::ParameterTree &config) {
    std::string key;
    if (config.hasKey("solver_type")) {
      key = "solver_type";
    } else if (config.hasSub("source_model")) {
      key = "source_model.type";
    }
    std::string feature_id = config.get<std::string>(key);
    if (get_experimental(feature_id) && (!enable_experimental_)) {
      config[key] = "experimental::" + config[key];
    }
    update_features(feature_id);
  }

  void update_features(const std::string feature_id) {
    relevant_features_.push_back(feature_id);
  }

  void print_citations() {
    relevant_features_.sort();
    relevant_features_.unique();
    relevant_features_.insert(relevant_features_.begin(), {"duneuro", "DUNE"});
    for (const std::string &feature_id : relevant_features_) {
      std::cout << get_description(feature_id) << "\n\n";
      std::vector<std::string> feature_papers = get_feature_papers(feature_id);
      for (const std::string &paper_id : feature_papers) {
        std::cout << get_citation(paper_id) << "\n\n";
      }
    }
  }

  bool get_experimental(const std::string &feature_id) const {
    return (experimental_features_.find(feature_id) !=
            experimental_features_.end());
  }

  std::vector<std::string>
  get_feature_papers(const std::string &feature_id) const {
    return features_.at(feature_id).feature_papers;
  }

  std::string get_description(const std::string &feature_id) const {
    return features_.at(feature_id).description;
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
