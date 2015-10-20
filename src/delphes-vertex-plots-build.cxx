#include <iostream>
#include <utility>
#include <vector>
#include <set>
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
#include "ndhist/Axis.hh"
#include "ndhist/Histogram.hh"

// === define constants ===
// const double GeV = 1;
// const double MeV = 0.001;

const double pi = std::atan2(0, -1);

// === list of histograms ===

const double VERTEX_MASS_MAX = 10;
const unsigned MAX_JETS = 20;
const unsigned MAX_TRACKS = 10;
const unsigned MAX_VERTEX = 10;

const float D0_RANGE = 1.0;
const float Z0_RANGE = 2.0;

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
  const Axis VXN  = {"vxn" , MAX_VERTEX + 1,  0.5, MAX_VERTEX + 1.5, ""};
  const Axis DRJV = {"drjv", 100, 0, 4};
  const Axis DPHI = {"dphi", 100, -pi, pi};
  const Axes vx_vars{
    PT, ETA, LXY, LSIG, EFRC, MASS, NTRK, DRJV, VXN, NVTX, DPHI};
  const Axes jet_vars{PT, ETA, NVTX};
  const Axes sv_vars{PT, ETA, LSIG, NVTX, NTRK, DRJV, MASS};

  const Axis TRK_WT = {"weight", 100, 0, 1};
  const Axis TRK_D0 = {"d0", 100, -D0_RANGE, D0_RANGE, "mm"};
  const Axis TRK_Z0 = {"z0", 100, -Z0_RANGE, Z0_RANGE, "mm"};
  const Axis TRK_D0SIG = {"d0sig", 100, -10, 10};
  const Axis TRK_Z0SIG = {"z0sig", 100, -10, 10};
  const Axis TRK_PT = {"pt", 100, 0, 100, "GeV"};
  const Axis TRK_DPHI = {"dphi", 100, -0.5, 0.5};
  const Axis TRK_DETA = {"deta", 100, -0.5, 0.5};
  const Axes trk_vars {
    TRK_WT, TRK_D0, TRK_Z0, TRK_D0SIG, TRK_Z0SIG, TRK_PT,
      TRK_DPHI, TRK_DETA};
  // const std::vector<Axis> all_vars{PT, ETA, LXY, EFRC, MASS, NTRK};
  const DMap jet_var_map(const Jet* jet);
  const DMap vx_var_map(const Jet& jet, const TSecondaryVertex& vxn);
  const DMap sv_var_map(const Jet& jet,
			  const std::vector<TSecondaryVertex>& vertices);
  const DMap trk_var_map(const TSecondaryVertexTrack& track);
}
std::ostream& operator<<(std::ostream& os, const var::DMap&);

// === misc utility functions ===
namespace {

  int fix_nan(var::DMap& map, const std::set<std::string>& silent) {
    int fixed = 0;
    for (auto& var: map) {
      if (isnan(var.second)) {
	if (silent.count(var.first)) {
	  fixed++;
	  var.second = -1;
	} else {
	  throw std::range_error(var.first + " is NaN");
	}
      }
    }
    return fixed;
  }

  std::map<std::string, std::vector<TSecondaryVertex> > sort_by_algo(
    const std::vector<TSecondaryVertex>& input) {
    std::map<std::string, std::vector<TSecondaryVertex> > out;
    for (auto vx: input) {
      std::string key = vx.config;
      key.erase(std::remove(key.begin(), key.end(), ':'), key.end());
      out[key].push_back(vx);
    }
    return out;
  }

  bool less_significant(const TSecondaryVertex& v1,
			const TSecondaryVertex& v2) {
    return v1.Lsig < v2.Lsig;
  }

  struct VxHists {
    AllPlanes vertex;
    AllPlanes tracks;
    AllPlanes tracks_hl;
    VxHists(var::Axes vxax, var::Axes trax):
      vertex(vxax), tracks(trax, true), tracks_hl(trax, true)
      {
      };
  };

  typedef std::map<int,std::map<std::string,std::map<int, AllPlanes>>> VxMap;
  typedef std::map<int,std::map<std::string, VxHists> > HlVarMap;
  typedef std::map<std::string, int> IMap;

  void fill_svx(VxMap& vx_map, HlVarMap& hl_map,
		const Jet& jet, IMap& counter) {
    // loop ovar algs
    for (const auto& algo_vx: sort_by_algo(jet.SecondaryVertices)){
      const auto& algo = algo_vx.first;
      auto vx_vec = algo_vx.second;
      // go from least to most significant
      std::sort(vx_vec.begin(), vx_vec.end(), less_significant);

      // loop over vertices
      const int n_vx = vx_vec.size();
      for (int iv = 0; iv < n_vx; iv++) {
	counter["n vtx"]++;
	auto vx_vars = var::vx_var_map(jet, vx_vec.at(iv));
	vx_vars[var::VXN.name] = iv + 1;
	vx_vars[var::NVTX.name] = n_vx;
	if (fix_nan(vx_vars, {var::LSIG.name})) counter["nan vx sig"]++;
	auto& vert_map = vx_map[jet.Flavor][algo];
	if (!vert_map.count(iv) ) {
	  vert_map.emplace(iv, AllPlanes(var::vx_vars, true));
	}
	auto& hists = vert_map.at(iv);
	hists.fill(vx_vars);
      } // end vx loop

      // fill high-level stuff
      auto sv_vars = var::sv_var_map(jet, vx_vec);
      fix_nan(sv_vars, {var::LSIG.name});
      // std::cout.width(40);
      // std::cout << algo << " flav: " << jet.Flavor << " " << sv_vars
      // 		<< std::endl;
      auto& flav_map = hl_map[jet.Flavor];
      if (!flav_map.count(algo)) {
	flav_map.emplace(algo, VxHists(var::sv_vars, var::trk_vars));
      }
      flav_map.at(algo).vertex.fill(sv_vars);

      // loop over constituent tracks
      auto& track_hists = flav_map.at(algo).tracks;
      for (const auto& vx: vx_vec) {
	for (const auto& trk: vx.tracks) {
	  track_hists.fill(var::trk_var_map(trk));
	}
      }
      // fill high-level tracks
      auto& hl_track_hists = flav_map.at(algo).tracks_hl;
      for (const auto& trk: jet.HLSecondaryVertexTracks) {
	hl_track_hists.fill(var::trk_var_map(trk));
      }
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
  // key by flavor, algo, vx number
  VxMap vx_map;
  // key by flavor, algo
  HlVarMap hl_map;
  std::map<std::string, int> counter;

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
      counter["total jets"]++;
      Jet* jet = root::as<Jet>(bJets->At(i_jet));
      auto const jet_vars = var::jet_var_map(jet);
      if (!jets_by_flavor.count(jet->Flavor)) {
	jets_by_flavor.emplace(jet->Flavor, var::jet_vars);
      }
      jets_by_flavor.at(jet->Flavor).fill(jet_vars);
      fill_svx(vx_map, hl_map, *jet, counter);
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
  auto by_flavor_algo = out_file.createGroup("algo");
  for (const auto& flav_itr: hl_map) {
    auto algo_group = by_flavor_algo.createGroup(
      std::to_string(flav_itr.first));
    for (const auto& algo_itr: flav_itr.second) {
      algo_itr.second.vertex.save_to(algo_group, algo_itr.first);
    }
  }
  auto tracks = out_file.createGroup("tracks");
  for (const auto& flav_itr: hl_map) {
    auto algo_group = tracks.createGroup(
      std::to_string(flav_itr.first));
    for (const auto& algo_itr: flav_itr.second) {
      algo_itr.second.tracks.save_to(algo_group, algo_itr.first);
      algo_itr.second.tracks_hl.save_to(algo_group, "high-level");
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

  const DMap trk_var_map(const TSecondaryVertexTrack& trk) {
    return {
      {TRK_WT.name, trk.weight},
      {TRK_D0.name, trk.d0},
      {TRK_Z0.name, trk.z0},
      {TRK_D0SIG.name, trk.d0 / trk.d0err},
      {TRK_Z0SIG.name, trk.z0 / trk.z0err},
      {TRK_PT.name, trk.pt},
      {TRK_DPHI.name, trk.dphi},
      {TRK_DETA.name, trk.deta}
    };
  }

  const DMap vx_var_map(const Jet& jet, const TSecondaryVertex& vx) {
    TVector3 jvec;
    jvec.SetPtEtaPhi(jet.PT, jet.Eta, jet.Phi);
    TVector3 vxvec(vx.x, vx.y, vx.z);
    double drjv = jvec.DeltaR(vxvec);
    double dphi = jvec.DeltaPhi(vxvec);
    return {
      {PT.name, jet.PT},
      {ETA.name, jet.Eta},
      {LXY.name, std::log1p(vx.Lxy)},
      {LSIG.name, std::log1p(vx.Lsig)},
      {EFRC.name, vx.eFrac},
      {MASS.name, vx.mass},
      {NTRK.name, vx.nTracks},
      {DRJV.name, drjv},
      {DPHI.name, dphi},
    };
  } // end vx_var_map
  const DMap sv_var_map(const Jet& jet,
			const std::vector<TSecondaryVertex>& vertices) {
    std::vector<TSecondaryVertex> over_sig_threshold;
    DMap output {
      {PT.name, jet.PT},
      {ETA.name, jet.Eta},
    };
    TVector3 jvec;
    jvec.SetPtEtaPhi(jet.PT, jet.Eta, jet.Phi);

    double sum_dr_tracks = 0;
    int sum_tracks = 0;
    int sum_vertices = 0;
    double sum_mass = 0;

    // copied these variables from jetfitter
    double sum_sig = 0;
    double sum_inverr = 0;
    const size_t n_vtx = vertices.size();
    for (size_t vxn = 1; vxn < n_vtx; vxn++) {
      const auto& vx = vertices.at(vxn);
      TVector3 vxv(vx.x, vx.y, vx.z);
      sum_vertices++;
      double delta_r = jvec.DeltaR(vxv);
      sum_dr_tracks += vx.nTracks * delta_r;
      sum_tracks += vx.nTracks;

      sum_mass += vx.mass;

      sum_sig += vxv.Mag() / vx.decayLengthVariance;
      sum_inverr += 1/vx.decayLengthVariance;
    }
    bool has_vx = (sum_vertices) > 0 && (sum_tracks > 0);
    output[LSIG.name] = has_vx ? std::log1p(sum_sig / sqrt(sum_inverr)) : 0;
    output[NVTX.name] = has_vx ? sum_vertices               : 0;
    output[NTRK.name] = has_vx ? sum_tracks                 : 0;
    output[DRJV.name] = has_vx ? sum_dr_tracks / sum_tracks : 3.0;
    output[MASS.name] = has_vx ? sum_mass                   : 0;
    return output;
  }

}

std::ostream& operator<<(std::ostream& os, const var::DMap& map) {
  using namespace std;
  for (const auto& item: map) {
    os << item.first << ": " << item.second << ", ";
  }
  return os;
}
