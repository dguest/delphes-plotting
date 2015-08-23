#include "root.hh"

#include "TFile.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include "misc_func.hh" 	// cli

#include "H5Cpp.h"
#include "ndhist/Axis.hh"
#include "ndhist/Histogram.hh"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>

const double pi = std::atan2(0, -1);

struct Hists
{
  Hists();
  void save(H5::CommonFG&);
  void save(H5::CommonFG&, const std::string&);
  void fill(const std::map<std::string, double>&);
  // tracking based
  Histogram transverse;
  Histogram longitudinal;
};

const size_t BINS = 200;
const unsigned flags = 0;

const Axis LONG   = {"z", BINS, -2, 2, "mm"};
const Axis TRANSX = {"x", BINS, -1, 1, "mm"};
const Axis TRANSY = {"y", BINS, -1, 1, "mm"};
const Axis TRANSR = {"r", BINS, -1, 1, "mm"};

Hists::Hists():
  transverse({TRANSX, TRANSY}),
  longitudinal({TRANSR, LONG})
{
}

void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.write_to(out_h5, #VAR)
  WRITE(transverse);
  WRITE(longitudinal);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}
void Hists::fill(const std::map<std::string, double>& vals) {
  transverse.fill(vals);
  longitudinal.fill(vals);
}

struct Point {
  double x;
  double y;
  double z;
  Point operator-(const Point& pt) {
    return {x - pt.x, y - pt.y, z - pt.z};
  }
  std::map<std::string, double> map() const {
    return {{"x", x}, {"y", y}, {"z", z}, {"r", std::hypot(x, y)}};
  }
};

Point alignWithJet(const Point& pt, const Jet& jet) {
  using namespace std;
  double jet_phi = jet.Phi;
  double r = hypot(pt.x, pt.y);
  double dphi = atan2(pt.y, pt.x) - jet_phi;
  Point out;
  if (r == 0) {
    out.x = 0;
    out.y = 0;
  } else {
    out.x = r*sin(dphi);
    out.y = r*cos(dphi);
  }
  out.z = pt.z;
  return out;
}

struct FlavorHists
{
  Hists bottom;
  Hists charm;
  Hists light;
  void fill(const Jet& jet);
  void save(H5::CommonFG& out);
};

namespace {
  void fill_one_reco_residuals(FlavorHists& hists, const Jet& jet);
}
void FlavorHists::fill(const Jet& jet) {
  fill_one_reco_residuals(*this, jet);
}
void FlavorHists::save(H5::CommonFG& out) {
  bottom.save(out, "bottom");
  charm.save(out, "charm");
  light.save(out, "light");
}

namespace {
  void fill_one_reco_residuals(FlavorHists& hists, const Jet& jet) {
    if (jet.SecondaryVertices.size() != 1) return;

    auto& sec = jet.SecondaryVertices.at(0);
    auto reco = alignWithJet({sec.X(), sec.Y(), sec.Z()}, jet);
    size_t n_vx = jet.TruthVertices.size();
    if (n_vx == 0) {
      hists.light.fill(reco.map());
    } else if (n_vx == 1 || n_vx == 2) {
      // TODO: is this the most significant, or second-most?
      auto& vx = jet.TruthVertices.at(0);
      auto tr = alignWithJet({vx.x, vx.y, vx.z}, jet);
      Point resid = reco - tr;
      if (n_vx == 1) hists.charm.fill(resid.map());
      else if (n_vx == 2) hists.bottom.fill(resid.map());
    }
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

  // Hists hists;
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

      // fill_hists(hists, *jet);
      flavhists.fill(*jet);
    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  flavhists.save(out_file);
  return 0;

}

