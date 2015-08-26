#include "root.hh"
#include "truth_tools.hh"

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

const Axis LONG   = {"z", BINS, -100, 100, "mm"};
const Axis TRANSX = {"x", BINS, -2, 2, "mm"};
const Axis TRANSY = {"y", BINS, -5, 10, "mm"};
const Axis TRANSR = {"r", BINS, 0, 20, "mm"};

Hists::Hists():
  transverse({TRANSX, TRANSY}),
  longitudinal({LONG, TRANSY})
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
template <class T>
Point as_point(const T& pt){
  return {pt.x, pt.y, pt.z};
}
template <>
Point as_point<TVector3>(const TVector3& pt){
  return {pt.X(), pt.Y(), pt.Z()};
}

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

struct ResponseHist
{
  ResponseHist();
  Histogram yproj;
  void fill(const Point& trth, const Point& reco);
  void save(H5::CommonFG& out, const std::string& subgrp);
};
ResponseHist::ResponseHist():
  yproj({{"true", BINS, -5, 15, "mm"}, {"reco", BINS, -5, 15, "mm"}})
{}
void ResponseHist::save(H5::CommonFG& out, const std::string& subgrp) {
  H5::Group grp(out.createGroup(subgrp));
  yproj.write_to(grp, "yproj");
}
void ResponseHist::fill(const Point& trth, const Point& reco) {
  yproj.fill({trth.y, reco.y});
}


template <class T>
struct FlavorHists
{
  T bottom;
  T charm;
  T light;
  void save(H5::CommonFG& out);
  void save(H5::CommonFG& out, std::string name);
};

template <class T>
void FlavorHists<T>::save(H5::CommonFG& out) {
  bottom.save(out, "bottom");
  charm.save(out, "charm");
  light.save(out, "light");
}
template <class T>
void FlavorHists<T>::save(H5::CommonFG& out, std::string name) {
  H5::Group out_group = out.createGroup(name);
  save(out_group);
}

namespace {
  void fill_one_reco_residuals(FlavorHists<Hists>& hists, const Jet& jet) {
    size_t nsv = jet.SecondaryVertices.size();
    if (nsv < 1) return;

    int max_trks_idx = 0;
    int max_tracks = 0;
    for (int vxn = 0; vxn < nsv; vxn++) {
      const auto& vx = jet.SecondaryVertices.at(vxn);
      if (vx.nTracks > max_tracks){
	max_tracks = vx.nTracks;
	max_trks_idx = vxn;
      }
    }
    const auto& sec = jet.SecondaryVertices.at(max_trks_idx);
    auto reco = alignWithJet({sec.X(), sec.Y(), sec.Z()}, jet);
    size_t n_vx = jet.TruthVertices.size();
    if (n_vx == 0) {
      hists.light.fill(reco.map());
    } else if (n_vx == 1 || n_vx == 2) {
      // vertices are ordered like the decay chain, first is first to decay
      auto& vx = jet.TruthVertices.at(0);
      auto tr = alignWithJet({vx.x, vx.y, vx.z}, jet);
      Point resid = reco - tr;
      auto map = resid.map();
      // for (auto& it: map) {
      // 	std::cout << it.first << " " << it.second << std::endl;
      // }
      if (n_vx == 1) hists.charm.fill(map);
      else if (n_vx == 2) hists.bottom.fill(map);
    }
  }
  void fill_truth(FlavorHists<Hists>& hists, const Jet& jet) {
    size_t n_vx = jet.TruthVertices.size();
    if (n_vx == 0) return;
    auto& vx = jet.TruthVertices.at(0);
    auto truth = alignWithJet(as_point(vx), jet);
    if (n_vx == 1) {
      hists.charm.fill(truth.map());
    } else if (n_vx == 2) {
      hists.bottom.fill(truth.map());
    }
  }
  void fill_reco(FlavorHists<Hists>& hists, const Jet& jet) {
    int flav = jet.Flavor;
    for (const auto& vx: jet.SecondaryVertices) {
      Point reco = alignWithJet(as_point<TVector3>(vx), jet);
      if (flav == 5) {
	hists.bottom.fill(reco.map());
      } else if (flav == 4) {
	hists.charm.fill(reco.map());
      } else {
	hists.light.fill(reco.map());
      }
    }
  }
  void fill_response(FlavorHists<ResponseHist>& hists, const Jet& jet) {
    if (jet.SecondaryVertices.size() == 0) return;
    for (const auto& vx: jet.TruthVertices) {
      int flav = truth::major_quark(vx.pdgid);
      Point reco = as_point<TVector3>(jet.SecondaryVertices.at(0));
      Point trth = as_point(vx);
      if (flav == 5) {
	hists.bottom.fill(trth, reco);
      } else if (flav == 4) {
	hists.charm.fill(trth, reco);
      } else {
	hists.light.fill(trth, reco);
      }
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
  FlavorHists<Hists> residuals;
  FlavorHists<Hists> truth;
  FlavorHists<Hists> reco;
  FlavorHists<ResponseHist> response;
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
      fill_one_reco_residuals(residuals,*jet);
      fill_truth(truth, *jet);
      fill_reco(reco, *jet);
      fill_response(response, *jet);
    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  residuals.save(out_file, "residuals");
  truth.save(out_file, "truth");
  reco.save(out_file, "reco");
  response.save(out_file, "response");

  return 0;

}

