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
#include "misc_func.hh"
#include "root.hh"
#include "truth_tools.hh"

#include "H5Cpp.h"
#include "Histogram.hh"

// === define constants ===
const double GeV = 1;
const double MeV = 0.001;

// b-tagging bits
namespace bit {
  // const unsigned B_TAG = 1 << 0;
  // const unsigned B_FLAVOR = 1 << 1;
}
const double pi = std::atan2(0, -1);

// === list of histograms ===

struct Hists {
  Hists();
  void save(std::string);
  void save(H5::CommonFG&);
  void save(H5::CommonFG&, const std::string&);
  Histogram n_tracks;
  Histogram n_jets;
  Histogram track_jet_dr;
  Histogram track_pt;
  Histogram track_d0;
  Histogram track_ip;
  Histogram particle_d0;
  Histogram particle_ip;
  Histogram track_z0;
  Histogram particle_z0;
  Histogram track_d0sig;
  Histogram track_z0sig;
  Histogram track_ipsig;
  Histogram initial_d0;
  Histogram track_d0phi;
  Histogram track_d0_dphi;
  Histogram jet_d0_dphi;
  Histogram smear_d0;
};

const unsigned MAX_TRACKS = 50;
const unsigned MAX_JETS = 20;
const float D0_RANGE = 1.0;
const float Z0_RANGE = 2.0;

Hists::Hists():
  n_tracks(MAX_TRACKS + 1, -0.5, MAX_TRACKS + 0.5),
  n_jets(MAX_JETS + 1, -0.5, MAX_JETS + 0.5),
  track_jet_dr(200, 0, 5),
  track_pt(200, 0, 200, "GeV"),
  track_d0(1000, -D0_RANGE, D0_RANGE, "mm"),
  track_ip(1000, -D0_RANGE, D0_RANGE, "mm"),
  particle_d0(1000, -D0_RANGE, D0_RANGE, "mm"),
  particle_ip(1000, -D0_RANGE, D0_RANGE, "mm"),
  track_z0(100, -Z0_RANGE, Z0_RANGE, "mm"),
  particle_z0(100, -Z0_RANGE, Z0_RANGE, "mm"),
  track_d0sig(1000, -30, 30, ""),
  track_z0sig(1000, -30, 30, ""),
  track_ipsig(1000, -30, 30, ""),
  initial_d0(1000, -D0_RANGE, D0_RANGE, "mm"),
  track_d0phi(1000, -pi, pi),
  track_d0_dphi(1000, -pi, pi),
  jet_d0_dphi(1000, -pi, pi),
  smear_d0(1000, -D0_RANGE, D0_RANGE, "mm")
{
}
void Hists::save(std::string output) {
  H5::H5File out_file(output, H5F_ACC_EXCL);
  save(out_file);
}
void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.write_to(out_h5, #VAR)
  WRITE(n_tracks);
  WRITE(n_jets);
  WRITE(track_jet_dr);
  WRITE(track_pt);
  WRITE(track_d0);
  WRITE(track_ip);
  WRITE(particle_d0);
  WRITE(particle_ip);
  WRITE(track_z0);
  WRITE(particle_z0);
  WRITE(track_d0sig);
  WRITE(track_z0sig);
  WRITE(track_ipsig);
  WRITE(initial_d0);
  WRITE(track_d0phi);
  WRITE(jet_d0_dphi);
  WRITE(track_d0_dphi);
  WRITE(smear_d0);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

// === misc utility functions ===

double get_ip(const Track* track, const Jet* jet) {
  // took this form the TrackCountingBTagging module
  // track perigees
  float xd = track->Xd;
  float yd = track->Yd;
  float d0 = track->Dxy;

  // jet momentum
  TLorentzVector jvec;
  jvec.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
  double jpx = jvec.Px();
  double jpy = jvec.Py();

  // signed impact parameter along jet axis
  double j_dot_d = jpx*xd + jpy*yd;
  // int sign = (j_dot_d > 0) ? 1 : -1;
  double ip = std::copysign(d0, j_dot_d);
  return ip;
}

void fill_track_hists(Hists& hists, const Track* track, const Jet* jet) {
  const Track* particle = root::as<Track>(track->Particle.GetObject());

  TLorentzVector jvec;
  jvec.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
  TLorentzVector trackvec;
  trackvec.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, 0);
  hists.track_jet_dr.fill(jvec.DeltaR(trackvec));
  hists.track_pt.fill(track->PT);

  double ip = get_ip(track, jet);
  hists.track_ip.fill(ip);
  double particle_ip = get_ip(particle, jet);
  // if (std::abs(particle_ip) > 1e-5) {
  hists.particle_ip.fill(particle_ip);
  // }

  float d0 = track->Dxy;
  hists.track_d0.fill(d0);
  float z0 = track->trkPar[TrackParam::Z0];
  hists.track_z0.fill(z0);
  {
    float d0_particle = particle->Dxy;
    // printf("d0 %f\n", d0_particle);
    hists.particle_d0.fill(d0_particle);
    float z0_particle = particle->trkPar[TrackParam::Z0];
    hists.particle_z0.fill(z0_particle);
    hists.smear_d0.fill(d0 - d0_particle);
  }
  float d0_cov = track->trkCov[TrackParam::D0D0];
  if (d0_cov > 0) {
    float d0sig = d0 / std::sqrt(d0_cov);
    hists.track_d0sig.fill(d0sig);
    float ipsig = ip / std::sqrt(d0_cov);
    hists.track_ipsig.fill(ipsig);
  }
  float z0_cov = track->trkCov[TrackParam::Z0Z0];
  if (z0_cov > 0) {
    float z0sig = z0 / std::sqrt(z0_cov);
    hists.track_z0sig.fill(z0sig);
  }
  GenParticle* gen = truth::get_gen_particle(track);
  if (!gen) {
    printf("none!\n");
    return;
  }
  double ix = gen->X;
  double iy = gen->Y;
  // std::cout << ix << " " << iy << " " << gen->PID << " "
  //        <<  " mothers " << gen->M1 << " " << gen->M2 <<std::endl;
  hists.initial_d0.fill(std::sqrt(ix*ix + iy*iy));
  hists.track_d0phi.fill(std::atan2(track->Y, track->X));
  {
    double xd = track->Xd;
    double yd = track->Yd;
    TLorentzVector d0vec(xd, yd, 0, 0);
    // printf("xd:  %f, yd: %f ", xd, yd);
    // printf("trx: %f, try: %f\n", trackvec.Px(), trackvec.Py());
    if (!std::isnan(xd) && !std::isnan(yd)) {
      hists.track_d0_dphi.fill(trackvec.DeltaPhi(d0vec));
      hists.jet_d0_dphi.fill(jvec.DeltaPhi(d0vec));
    }
  }
}

// === command line interface ===

struct CLI
{
  std::string out_name;
  std::string input;
  int err_code;
  void usage(std::string prname) {
    std::cerr << "usage: " << prname << ": <input> [<output>]" << std::endl;
  }
  int read(int argc, char* argv[]) {
    using namespace std;
    if (argc == 1 || argc > 3) {
      usage(argv[0]);
      err_code = 1;
      return err_code;
    }

    input = argv[1];
    if (!exists(input)) {
      cerr << input << " doesn't exist, exiting!" << endl;
      err_code = 1;
      return 1;
    }

    if (argc > 2) {
      out_name = argv[2];
    } else {
      out_name = "test.h5";
    }
    if (exists(out_name)) {
      cerr << out_name << " exists, exiting" << endl;
      err_code = 1;
      return err_code;
    }
    return 0;
  }
};

// ____________________________________________________
// main function


int main(int argc, char *argv[])
{
  gROOT->SetBatch();

  CLI cli;
  cli.read(argc, argv);
  if (cli.err_code) return cli.err_code;

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(cli.input.c_str());

  // Create object of class ExRootTreeReader
  ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
  long long int numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  treeReader->UseBranch("Track");
  TClonesArray* bPart = treeReader->UseBranch("Particle");
  treeReader->UseBranch("OriginalTrack");
  treeReader->UseBranch("Electron");
  treeReader->UseBranch("Photon");
  treeReader->UseBranch("Muon");
  treeReader->UseBranch("Tower");

  treeReader->UseBranch("EFlowTrack");
  // treeReader->UseBranch("EFlowPhoton");
  // treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray* bJets = treeReader->UseBranch("Jet");

  Hists hists;
  Hists b_jet_hists;
  Hists light_jet_hists;
  Hists leading_track_b_hists;
  Hists leading_track_light_hists;
  Hists matched_track_hists;
  Hists high_pt_tracks;
  Hists low_pt_tracks;
  std::map<int, int> b_decay_pids;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int n_jets = bJets->GetEntries();
    hists.n_jets.fill(n_jets);

    // loop over jets
    for (int i_jet = 0; i_jet < n_jets; i_jet++) {
      Jet* jet = root::as<Jet>(bJets->At(i_jet));
      int n_constituents = jet->Constituents.GetEntriesFast();

      // bool b_label = (jet->BTag & bit::B_FLAVOR);
      bool b_label = jet->Flavor == 5;

      // loop over all the tracks
      int n_tracks = 0;
      std::map<double, Track*> track_by_pt;
      // if (b_label) puts("new b-jet");
      // else puts("light jet");
      for (int iii = 0; iii < n_constituents; iii++) {
        TObject* obj = jet->Constituents.At(iii);
        if (!root::is<Track>(obj)) continue;
        auto* track = root::as<Track>(obj);
        n_tracks++;
        track_by_pt[track->PT] = track;

        // walk truth record, see if this track comes from something
        // interesting
        std::vector<std::deque<int> > sequences = {
          {5, 6}, {-5, -6},
          {5, -6}, {-5, 6},     // how does this happen?
          {5, 25}, {-5, 25} };
        std::vector<GenParticle*> daughters;
        for (const auto& seq: sequences) {
          auto daut = truth::get_parent(track, bPart, seq, 2);
          if (daut) daughters.push_back(daut);
        }
        if (daughters.size() == 1){
          const auto* daut = daughters.back();
          using std::pow;
          using std::sqrt;
          b_decay_pids[daut->PID]++;
          fill_track_hists(matched_track_hists, track, jet);
        } else if (daughters.size() > 1) {
          Track* particle = root::as<Track>(track->Particle.GetObject());
          std::cout << "found multiple matches for particle "
                    << "part d0: " << particle->Dxy << " "
                    << particle->trkPar[TrackParam::D0] << " "
                    << "trk d0: " << track->Dxy << " "
                    << track->trkPar[TrackParam::D0] << " "
                    << std::endl;
          for (const auto* daut: daughters) {
          // === dump some debugging info ===
          std::cout << "part " << truth::map_particle(daut->PID) << " "
                    << daut->PID << " "
                    << sqrt(pow(daut->X,2) + pow(daut->Y,2)) << " "
                    << std::endl;
          }
        }

        // fill hists
        fill_track_hists(hists, track, jet);
        if (b_label) {
          fill_track_hists(b_jet_hists, track, jet);
        } else {
          fill_track_hists(light_jet_hists, track, jet);
        }

	// low and high pt tracks for comp to
	// ATLAS-CONF-2014-047
	double track_theta = track->trkPar[TrackParam::THETA];
	double trk_eff_pt = track->PT * std::sqrt(track_theta);
	if (400*MeV < trk_eff_pt && trk_eff_pt < 500*MeV) {
	  fill_track_hists(low_pt_tracks, track, jet);
	} else if (track->PT > 20*GeV) {
	  fill_track_hists(high_pt_tracks, track, jet);
	}
      } // end loop over jet tracks
      hists.n_tracks.fill(n_tracks);
      if (track_by_pt.size() > 0) {
        const auto* track = track_by_pt.crbegin()->second;
        if (b_label) fill_track_hists(leading_track_b_hists, track, jet);
        else fill_track_hists(leading_track_light_hists, track, jet);
      }
    } // end loop over jets
  }   // end loop over events
  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  hists.save(out_file, "all_jets");
  b_jet_hists.save(out_file, "b_jets");
  light_jet_hists.save(out_file,"light_jets");
  leading_track_light_hists.save(out_file, "leading_track_light");
  leading_track_b_hists.save(out_file, "leading_track_b");
  matched_track_hists.save(out_file, "matched_track_hists");
  high_pt_tracks.save(out_file, "high_pt_tracks");
  low_pt_tracks.save(out_file, "low_pt_tracks");

  // dump list of most common particle decays
  std::vector<std::pair<int,int>> number_and_pid;
  for (const auto& itr: b_decay_pids) {
    number_and_pid.emplace_back(itr.second,itr.first);
  }
  std::sort(number_and_pid.begin(), number_and_pid.end());
  for (const auto& itr: number_and_pid) {
    printf("%s, num: %i\n",
           truth::map_particle(itr.second).c_str(), itr.first);
  }
  return 0;

}
