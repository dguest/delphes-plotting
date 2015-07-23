#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <deque>

#include "root.hh"
#include "AllPlanes.hh"

#include "TFile.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include "misc_func.hh" 	// cli
// #include "truth_tools.hh"

#include "H5Cpp.h"
#include "Axis.hh"
#include "Histogram.hh"

// === define constants ===
// const double GeV = 1;
// const double MeV = 0.001;

const double pi = std::atan2(0, -1);

// === list of histograms ===

const double VERTEX_MASS_MAX = 10;
const unsigned MAX_JETS = 20;
const unsigned MAX_TRACKS = 30;
const unsigned MAX_VERTEX = 10;

// all the 2d hist stuff
namespace var {
  typedef std::vector<Axis> Axes;
  typedef std::map<std::string, double> DMap;
  // names for variables
  const Axis PT   = {"pt"  , 100, 0, 100, "Gev"};
  const Axis ETA  = {"eta" , 100, -2.7, 2.7, ""};
  const Axis LXY  = {"lxy" , 100, 0, 7, "log1p mm"};
  const Axis LSIG = {"lsig", 100, 0, 7, "log1p"};
  const Axis EFRC = {"efrc", 100, 0, 1.001, ""};
  const Axis MASS = {"mass", 100, 0, 15, "GeV"};
  const Axis NTRK = {"ntrk", MAX_TRACKS + 1, -0.5, MAX_TRACKS + 0.5, ""};
  const Axis NVTX = {"nvtx", MAX_VERTEX + 1, -0.5, MAX_VERTEX + 0.5, ""};
  const Axis VXN  = {"vxn" , MAX_VERTEX + 1, -0.5, MAX_VERTEX + 0.5, ""};
  const Axis DRJV = {"drjv", 100, 0, 4};
  const Axes vx_vars{PT, ETA, LXY, LSIG, EFRC, MASS, NTRK, DRJV, VXN};
  const Axes jet_vars{PT, ETA, NVTX};
  // const std::vector<Axis> all_vars{PT, ETA, LXY, EFRC, MASS, NTRK};
  const DMap jet_var_map(const Jet* jet);
  const DMap vx_var_map(const Jet& jet, const SecondaryVertex& vxn);
}

// === misc utility functions ===
namespace {
  std::map<std::string, std::vector<SecondaryVertex> > sort_by_algo(
    const std::vector<SecondaryVertex>& input) {
    std::map<std::string, std::vector<SecondaryVertex> > out;
    for (auto vx: input) {
      out[vx.config].push_back(vx);
    }
    return out;
  }

  bool less_significant(const SecondaryVertex& v1, const SecondaryVertex& v2) {
    return v1.Lsig < v2.Lsig;
  }

  typedef std::map<int,std::map<std::string,std::map<int, AllPlanes>>> VxMap;
  typedef std::map<std::string, int> IMap;

  void fill_svx(VxMap& vx_map, const Jet& jet, IMap& counter) {
    for (const auto& algo_vx: sort_by_algo(jet.SecondaryVertices)){
      const auto& algo = algo_vx.first;
      auto vx_vec = algo_vx.second;
      std::sort(vx_vec.begin(), vx_vec.end(), less_significant);
      // go from most to least significant
      std::reverse(vx_vec.begin(), vx_vec.end());
      const int n_vx = vx_vec.size();
      for (int iv = 0; iv < n_vx; iv++) {
	counter["n vtx"]++;
	auto vx_vars = var::vx_var_map(jet, vx_vec.at(iv));
	vx_vars[var::VXN.name] = iv;
	if (isnan(vx_vars[var::LSIG.name])) {
	  counter["nan vx sig"]++;
	  vx_vars[var::LSIG.name] = -1;
	}
	auto& vert_map = vx_map[jet.Flavor][algo];
	if (!vert_map.count(iv) ) {
	  vert_map.emplace(iv, var::vx_vars);
	}
	auto& hists = vert_map.at(iv);
	hists.fill(vx_vars);
      } // end vx loop
    }   // end algo loop
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

  using namespace std;
  std::map<int, AllPlanes> jets_by_flavor;
  // key be flavor, algo, vx number
  VxMap vx_map;
  std::map<std::string, int> counter;

  // Loop over all events
  std::cout << "looping over " << numberOfEntries << " entries" << std::endl;
  int onem = numberOfEntries / 1000;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    if (entry % onem == 0) std::cout << entry << " processed\r" << std::flush;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int n_jets = bJets->GetEntries();

    // loop over jets
    for (int i_jet = 0; i_jet < n_jets; i_jet++) {
      counter["total jets"]++;
      Jet* jet = root::as<Jet>(bJets->At(i_jet));
      auto const jet_vars = var::jet_var_map(jet);
      if (!jets_by_flavor.count(jet->Flavor)) {
	jets_by_flavor.emplace(jet->Flavor, var::jet_vars);
      }
      jets_by_flavor.at(jet->Flavor).fill(jet_vars);
      fill_svx(vx_map, *jet, counter);
    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  auto by_flavor = out_file.createGroup("jets");
  for (const auto& fl_pls: jets_by_flavor){
    fl_pls.second.save_to(by_flavor, std::to_string(fl_pls.first));
  }
  auto by_flavor_vxn = out_file.createGroup("vx");
  for (const auto& flav_itr: vx_map) {
    auto algo_group = by_flavor_vxn.createGroup(
      std::to_string(flav_itr.first));
    for (const auto& algo_itr: flav_itr.second) {
      auto vx_group = algo_group.createGroup(algo_itr.first);
      for (const auto& vx_itr: algo_itr.second) {
	vx_itr.second.save_to(vx_group, std::to_string(vx_itr.first));
      }
    }
  }

  for (const auto& st_ct: counter) {
    std::cout << st_ct.first << ": " << st_ct.second
	      << std::endl;
  }
  return 0;

}

namespace var {
  const DMap jet_var_map(const Jet* jet) {
    auto sorted = sort_by_algo(jet->SecondaryVertices);
    size_t max_vx = 0;
    for (auto alg: sorted) {
      max_vx = std::max(alg.second.size(), max_vx);
    }
    return {
      {PT.name, jet->PT},
      {ETA.name, jet->Eta},
      {NVTX.name , max_vx},
    };
  } // end jet_var_map
  const DMap vx_var_map(const Jet& jet, const SecondaryVertex& vx) {
    TVector3 jvec;
    jvec.SetPtEtaPhi(jet.PT, jet.Eta, jet.Phi);
    double drjv = jvec.DeltaR(vx);
    return {
      {PT.name, jet.PT},
      {ETA.name, jet.Eta},
      {LXY.name, std::log1p(vx.Lxy)},
      {LSIG.name, std::log1p(vx.Lsig)},
      {EFRC.name, vx.eFrac},
      {MASS.name, vx.mass},
      {NTRK.name, vx.nTracks},
      {DRJV.name, drjv}
    };
  } // end vx_var_map
}
