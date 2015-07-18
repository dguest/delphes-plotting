#include "AllPlanes.hh"

#include "Histogram.hh"

#include "H5Cpp.h"

#include <utility>

AllPlanes::AllPlanes(const std::vector<Axis>& list)
{
  int n_ax = list.size();
  for (int iii = 0; iii < n_ax; iii++) {
    for (int jjj = iii + 1; jjj < n_ax; jjj++) {
      const auto& ax1 = list.at(iii);
      const auto& ax2 = list.at(jjj);
      std::string full_name = ax1.name + "_" + ax2.name;
      // Histogram hist({ax1, ax2} ,hist::flat_attributes);
      Histogram hist({ax1, ax2} );
      _hists.push_back(std::make_pair(full_name, hist));
    }
  }
}
void AllPlanes::fill(const std::map<std::string, double>& variables) {
  for (auto& nhist: _hists) {
    nhist.second.fill(variables);
  }
}
void AllPlanes::save_to(H5::CommonFG& group, std::string subgroup) {
  H5::Group out_group = group.createGroup(subgroup);
  save_to(out_group);
}
void AllPlanes::save_to(H5::CommonFG& group) {
  for (const auto& nhist: _hists) {
    nhist.second.write_to(group, nhist.first);
  }
}
