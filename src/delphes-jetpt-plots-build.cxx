#include <iostream>
// #include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
// #include <deque>

#include "root.hh"
// #include "AllPlanes.hh"

#include "TFile.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include "misc_func.hh" 	// cli
// #include "truth_tools.hh"

#include "H5Cpp.h"
#include "ndhist/Axis.hh"
#include "ndhist/Histogram.hh"

// === define constants ===
// const double GeV = 1;
// const double MeV = 0.001;

const double pi = std::atan2(0, -1);

// === list of histograms ===

// const double VERTEX_MASS_MAX = 10;
// const unsigned MAX_JETS = 20;
// const unsigned MAX_TRACKS = 30;
// const unsigned MAX_VERTEX = 10;

// === misc utility functions ===

struct Hists
{
  Hists();
  void save(H5::CommonFG&);
  void save(H5::CommonFG&, const std::string&);
  // tracking based
  Histogram jet_pt;
};

// const float D0_RANGE = 1.0;
// const float Z0_RANGE = 2.0;
const size_t BINS = 500;
const unsigned flags = 0;

Hists::Hists():
  jet_pt({{"jet_pt", BINS, 0, 1000, "GeV"}}, flags)
{
}

void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.write_to(out_h5, #VAR)
  WRITE(jet_pt);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

namespace {
  void fill_hists(Hists& hists, const Jet& jet) {
    hists.jet_pt.fill(jet.PT);
  }
}

struct FlavorHists
{
  Hists bottom;
  Hists charm;
  Hists light;
  void fill(const Jet& jet);
  void save(H5::CommonFG& out);
  void save(H5::CommonFG& out, std::string subdir);
};

void FlavorHists::fill(const Jet& jet) {
  unsigned ftl = jet.Flavor;
  if (ftl == 5) fill_hists(bottom, jet);
  else if (ftl == 4) fill_hists(charm, jet);
  else if (ftl < 4 || ftl == 21) fill_hists(light, jet);
  else assert(false);
}
void FlavorHists::save(H5::CommonFG& out) {
  bottom.save(out, "bottom");
  charm.save(out, "charm");
  light.save(out, "light");
}
void FlavorHists::save(H5::CommonFG& out, std::string subdir) {
  auto subgroup = out.createGroup(subdir);
  save(subgroup);
}

// ____________________________________________________
// main function


int main(int argc, char *argv[])
{
  gROOT->SetBatch();

  FileCLI cli(argc, argv);

  // Create chain of root trees
  TChain chain("Delphes");
  for (const auto& in: cli.in_files()) {
    chain.Add(in.c_str());
  }

  // Create object of class ExRootTreeReader
  ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
  long long int numberOfEntries = treeReader->GetEntries();
  TClonesArray* bJets = treeReader->UseBranch("Jet");

  FlavorHists flavhists;
  using namespace std;

  // Loop over all events
  std::cout << "looping over " << numberOfEntries << " entries" << std::endl;
  int onem = std::max<int>(numberOfEntries / 1000, 1);
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    if (entry % onem == 0) std::cout << entry << " processed\r" << std::flush;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int n_jets = bJets->GetEntries();

    // loop over jets
    for (int i_jet = 0; i_jet < n_jets; i_jet++) {
      Jet* jet = root::as<Jet>(bJets->At(i_jet));
      if (std::abs(jet->Eta) > 2.5) continue;
      if (std::abs(jet->PT) < 20) continue;

      flavhists.fill(*jet);
    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_file(), H5F_ACC_EXCL);
  flavhists.save(out_file, "jet_vars");
  return 0;

}

