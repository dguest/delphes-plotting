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
  Histogram jetProb;
  Histogram track2d0;
  Histogram track3d0;

  // vertex based
  Histogram lsig;
  Histogram drjet;
  Histogram mass;
};

// const float D0_RANGE = 1.0;
// const float Z0_RANGE = 2.0;
const size_t BINS = 1000;
const unsigned flags = 0;

Hists::Hists():
  jetProb(BINS, -20, 0, "log", flags),
  track2d0(BINS, -10, 10, "", flags),
  track3d0(BINS, -10, 10, "", flags),
  lsig(BINS, 0, 7, "log1p", flags),
  drjet(BINS, 0, 4, "", flags),
  mass(BINS, 0, 50, "GeV", flags)
{
}

void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.write_to(out_h5, #VAR)
  WRITE(jetProb);
  WRITE(track2d0);
  WRITE(track3d0);
  WRITE(lsig);
  WRITE(drjet);
  WRITE(mass);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

namespace {
  void fill_hists(Hists& hists, const Jet& jet) {
    hists.jetProb.fill(jet.jetProb >= 0 ? std::log(jet.jetProb): -1);
    hists.track2d0.fill(jet.track2d0sig);
    hists.track3d0.fill(jet.track3d0sig);
    hists.lsig.fill(jet.svLsig > 0 ? std::log1p(jet.svLsig): -1);
    hists.drjet.fill(jet.svDrJet);
    hists.mass.fill(jet.svMass);
  }
}
// ____________________________________________________
// main function


int main(int argc, char *argv[])
{
  gROOT->SetBatch();

  CLI cli;
  if (cli.read(argc, argv)) return cli.err_code();

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(cli.in_name.c_str());

  // Create object of class ExRootTreeReader
  ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
  long long int numberOfEntries = treeReader->GetEntries();
  TClonesArray* bJets = treeReader->UseBranch("Jet");

  Hists hists;
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

      fill_hists(hists, *jet);
    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  hists.save(out_file, "high-level");
  return 0;

}

