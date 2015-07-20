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
  void save_to(H5::CommonFG& , std::string subgroup) const;
  void save_to(H5::CommonFG&) const;
private:
  std::vector<std::pair<std::string,Histogram> > _hists;
};

#endif
