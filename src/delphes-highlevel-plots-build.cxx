#include "parse_json.hh"
#include "misc_func.hh" 	// cli
#include "LightweightNeuralNetwork.hh"
#include "NanReplacer.hh"
#include "hl_var_map.hh"
// #include "truth_tools.hh"
#include "root.hh"

#include "TFile.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include "H5Cpp.h"
#include "ndhist/Axis.hh"
#include "ndhist/Histogram.hh"

#include <iostream>
// #include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <fstream>
// #include <deque>

// #include "AllPlanes.hh"


// === define constants ===
// const double GeV = 1;
// const double MeV = 0.001;

const double pi = std::atan2(0, -1);
static const char* JSON_NN_FILE = "pheno-wt.json";
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
  Histogram track2z0;
  Histogram track3z0;
  Histogram n_sig_track;
  Histogram width_eta;
  Histogram width_phi;

  // vertex based
  Histogram lsig;
  Histogram n_vx_track;
  Histogram drjet;
  Histogram mass;
  Histogram nsecvtx;
  Histogram efrac;

  lwt::LightweightNeuralNetwork* lwtnn;
  Histogram julian;
  lwt::NanReplacer* nan_rep;
};

// const float D0_RANGE = 1.0;
// const float Z0_RANGE = 2.0;
const size_t BINS = 10000;
const unsigned flags = hist::eat_nan;
const unsigned MAX_VERTEX = 4;
const unsigned MAX_TRACKS = 25;

Hists::Hists():
  jetProb({{"jetProb", BINS, -50, 0, "log"}}, flags),
  track2d0({{"track2d0sig", BINS, -2, 20, ""}}, flags),
  track3d0({{"track3d0sig", BINS, -2, 20, ""}}, flags),
  track2z0({{"track2z0sig", BINS, -2, 20, ""}}, flags),
  track3z0({{"track3z0sig", BINS, -2, 20, ""}}, flags),
  n_sig_track({{"n_sig_track", MAX_TRACKS + 1, -0.5, MAX_TRACKS + 0.5, ""}},
	      flags),
  width_eta({{"width_eta", BINS, 0, 0.5, ""}}, flags),
  width_phi({{"width_phi", BINS, 0, 0.5, ""}}, flags),

  lsig({{"lsig", BINS, 0, 7, "log1p"}}, flags),
  n_vx_track({{"n_vx_track", MAX_TRACKS + 1, -0.5, MAX_TRACKS + 0.5, ""}},
	     flags),
  drjet({{"drjet", BINS, 0, 10, ""}}, flags),
  mass({{"mass", BINS, 0, 10, "GeV"}}, flags),
  nsecvtx({{"nsecvtx", MAX_VERTEX + 1, -0.5, MAX_VERTEX + 0.5, ""}}, flags),
  efrac({{"efrc", BINS, 0, 1.00001, ""}}, flags),
  lwtnn(0),
  julian({{"discriminant", BINS, 0, 1.00001, ""}}, flags),
  nan_rep(0)
{
}

void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.write_to(out_h5, #VAR)
  WRITE(jetProb);
  WRITE(track2d0);
  WRITE(track3d0);
  WRITE(track2z0);
  WRITE(track3z0);
  WRITE(n_sig_track);
  WRITE(width_eta);
  WRITE(width_phi);

  WRITE(lsig);
  WRITE(n_vx_track);
  WRITE(drjet);
  WRITE(mass);
  WRITE(nsecvtx);
  WRITE(efrac);

  if (lwtnn) WRITE(julian);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

namespace {
  void fill_hists(Hists& hists, const Jet& jet, bool do_ml) {
    hists.jetProb.fill(jet.jetProb >= 0 ? std::log(jet.jetProb): -1);
    hists.track2d0.fill(jet.track2d0sig);
    hists.track3d0.fill(jet.track3d0sig);
    hists.track2z0.fill(jet.track2z0sig);
    hists.track3z0.fill(jet.track3z0sig);
    hists.n_sig_track.fill(jet.tracksOverIpThreshold);
    hists.width_eta.fill(jet.jetWidthEta);
    hists.width_phi.fill(jet.jetWidthPhi);

    const auto& hl_svx = do_ml ?
      jet.MLSecondaryVertex : jet.HLSecondaryVertex;
    hists.lsig.fill(hl_svx.svLsig > 0 ? std::log1p(hl_svx.svLsig): -1);
    hists.n_vx_track.fill(hl_svx.svNTracks < 0 ? 0 : hl_svx.svNTracks);
    hists.drjet.fill(hl_svx.svDrJet);
    hists.mass.fill(hl_svx.svMass);
    hists.nsecvtx.fill(hl_svx.svNVertex < 0 ? 0 : hl_svx.svNVertex);
    hists.efrac.fill(hl_svx.svEnergyFraction);

    if (hists.lwtnn) {
      const auto hl_vars = hists.nan_rep->replace(hl_var_map(jet));
      // for (const auto& var: hl_vars) {
      //   std::cout << var.first << " " << var.second << " ";
      // }
      // std::cout << "\n";
      const auto disc = hists.lwtnn->compute(hl_vars);
      hists.julian.fill(disc.at("out_0"));
    }
  }
}

struct FlavorHists
{
  Hists bottom;
  Hists charm;
  Hists light;
  void fill(const Jet& jet, bool do_ml);
  void save(H5::CommonFG& out);
  void save(H5::CommonFG& out, std::string subdir);
  void set_nn(lwt::LightweightNeuralNetwork* network,
              lwt::NanReplacer* replacer);
};

void FlavorHists::fill(const Jet& jet, bool do_ml) {
  unsigned ftl = jet.Flavor;
  if (ftl == 5) fill_hists(bottom, jet, do_ml);
  else if (ftl == 4) fill_hists(charm, jet, do_ml);
  else if (ftl < 4 || ftl == 21) fill_hists(light, jet, do_ml);
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
void FlavorHists::set_nn(lwt::LightweightNeuralNetwork* network,
                         lwt::NanReplacer* rep) {
  bottom.lwtnn = network;
  charm.lwtnn = network;
  light.lwtnn = network;
  bottom.nan_rep = rep;
  charm.nan_rep =  rep;
  light.nan_rep =  rep;
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

  Hists hists;
  FlavorHists flavhists;
  FlavorHists flavhists_ml;
  using namespace std;

  if (exists(JSON_NN_FILE)) {
    std::cout << "loading nn file" << std::endl;
    std::ifstream file(JSON_NN_FILE);
    lwt::JSONConfig config = lwt::parse_json(file);
    auto* hl_nn = new lwt::LightweightNeuralNetwork(config.inputs,
                                                    config.layers,
                                                    config.outputs);
    auto* rep = new lwt::NanReplacer(config.defaults, lwt::rep::all);
    flavhists.set_nn(hl_nn, rep);
  }

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

      fill_hists(hists, *jet, false);
      flavhists.fill(*jet, false);
      flavhists_ml.fill(*jet, true);
    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_file(), H5F_ACC_TRUNC);
  hists.save(out_file, "all");
  flavhists.save(out_file, "high-level");
  flavhists_ml.save(out_file, "med-level");
  return 0;

}

