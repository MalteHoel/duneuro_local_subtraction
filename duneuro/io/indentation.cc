#include <config.h>
#include <duneuro/io/indentation.hh>

// implementation of indentation class
namespace duneuro {
  // class for easier handling of indentations

  Indentation::Indentation(size_t starting_level)
  {
    indentation_level_ = starting_level;
  }

  Indentation& Indentation::operator++() 
  {
    ++indentation_level_;
    return *this;
  }
  
  Indentation& Indentation::operator--()
  {
    --indentation_level_;
    return *this;
  }
  
  Indentation Indentation::operator++(int)
  {
    Indentation old(indentation_level_);
    ++indentation_level_;
    return old;
  }
  
  Indentation Indentation::operator--(int)
  {
    Indentation old(indentation_level_);
    --indentation_level_;
    return old;
  }

  std::ostream& operator<<(std::ostream& out, const Indentation& indentation)
  {
    for(size_t i = 0; i < indentation.indentation_level_; ++i) {
      out << "\t";
    }
    return out;
  }
  
  std::ofstream& operator<<(std::ofstream& out, const Indentation& indentation)
  {
    for(size_t i = 0; i < indentation.indentation_level_; ++i) {
      out << "\t";
    }
    
    return out;
  }
} // end namespace duneuro
