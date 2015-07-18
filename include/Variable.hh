#ifndef VARIABLE_HH
#define VARIABLE_HH

#include <string>

struct Variable {
  std::string name;
  int bins;
  double min;
  double max;
  std::string units;
};

#endif
