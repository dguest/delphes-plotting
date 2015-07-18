#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <deque>

#include "TROOT.h"
#include "TSystem.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"

#include "ExRootTreeReader.h"
#include "misc_func.hh" 	// cli
#include "root.hh"
#include "truth_tools.hh"

#include "H5Cpp.h"
#include "Histogram.hh"

// === define constants ===
// const double GeV = 1;
// const double MeV = 0.001;

const double pi = std::atan2(0, -1);

// === list of histograms ===

struct Hists {
  Hists();
  void save(std::string);
  void save(H5::CommonFG&);
  void save(H5::CommonFG&, const std::string&);
  Histogram vertex_mass;
  Histogram n_jets;
};

const double VERTEX_MASS_MAX = 10;
const unsigned MAX_JETS = 20;

Hists::Hists():
  vertex_mass(100, 0, VERTEX_MASS_MAX, "GeV"),
  n_jets(MAX_JETS + 1, -0.5, MAX_JETS + 0.5)
{
}
void Hists::save(std::string output) {
  H5::H5File out_file(output, H5F_ACC_EXCL);
  save(out_file);
}
void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.write_to(out_h5, #VAR)
  WRITE(vertex_mass);
  WRITE(n_jets);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

// === misc utility functions ===

void fill_vx_hists(Hists& hists, const Jet* jet) {
  TLorentzVector jvec;
  jvec.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
  hists.vertex_mass.fill(jet->SecondaryVertexMass);
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

  // Loop over all events
  std::cout << "looping over " << numberOfEntries << " entries" << std::endl;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int n_jets = bJets->GetEntries();
    hists.n_jets.fill(n_jets);

    // loop over jets
    for (int i_jet = 0; i_jet < n_jets; i_jet++) {
      Jet* jet = root::as<Jet>(bJets->At(i_jet));
      bool b_label = jet->Flavor == 5;
      if (b_label) {
	fill_vx_hists(hists, jet);
      }
    } // end loop over jets
  }   // end loop over events
  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  hists.save(out_file, "jet_vx");

  return 0;

}
