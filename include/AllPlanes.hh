#ifndef ALL_PLANES_HH
#define ALL_PLANES_HH

#include "Axis.hh"
#include "Histogram.hh"

#include <vector>
#include <map>
#include <utility>

namespace H5 {
  class CommonFG;
}

class AllPlanes
{
public:
  AllPlanes(const std::vector<Axis>& list);
  void fill(const std::map<std::string, double>& variables);
  void save_to(H5::CommonFG& , std::string subgroup);
  void save_to(H5::CommonFG&);
private:
  std::vector<std::pair<std::string,Histogram> > _hists;
};

#endif
