#ifndef DUNEURO_EEG_SOURCE_MODELS_INDENTATION_HH
#define DUNEURO_EEG_SOURCE_MODELS_INDENTATION_HH

#include <iostream>                                   // include for writing to std::cout
#include <fstream>                                    // include for writing to files 

namespace duneuro {
  // class for easier handling of indentations
  class Indentation
  {
  public:
    explicit Indentation(size_t starting_level = 0);
    Indentation& operator++();
    Indentation& operator--();
    Indentation operator++(int);
    Indentation operator--(int);
    friend std::ostream& operator<<(std::ostream& out, const Indentation& indentation);
    friend std::ofstream& operator<<(std::ofstream& out, const Indentation& indentation);
  private:
    size_t indentation_level_;
  };
  
  std::ostream& operator<<(std::ostream& out, const Indentation& indentation);
  std::ofstream& operator<<(std::ofstream& out, const Indentation& indentation);
}
#endif // DUNEURO_EEG_SOURCE_MODELS_INDENTATION_HH
