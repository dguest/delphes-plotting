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

const double GeV = 1;
const double MeV = 0.001;

// b-tagging bits
namespace bit {
  const unsigned B_TAG = 1 << 0;
  const unsigned B_FLAVOR = 1 << 1;
}

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
};

const unsigned MAX_TRACKS = 200;
const unsigned MAX_JETS = 20;
const float D0RNG = 0.1;
const float Z0RNG = 0.1;

Hists::Hists():
  n_tracks(MAX_TRACKS, -0.5, MAX_TRACKS + 0.5),
  n_jets(MAX_JETS, -0.5, MAX_JETS + 0.5),
  track_jet_dr(200, 0, 5),
  track_pt(200, 0, 200, "GeV"),
  track_d0(100, -D0RNG, D0RNG, "mm"),
  track_ip(100, -D0RNG, D0RNG, "mm"),
  particle_d0(100, -D0RNG, D0RNG, "mm"),
  particle_ip(100, -D0RNG, D0RNG, "mm"),
  track_z0(100, -Z0RNG, Z0RNG, "mm"),
  particle_z0(100, -Z0RNG, Z0RNG, "mm"),
  track_d0sig(1000, -30, 30, ""),
  track_z0sig(1000, -30, 30, ""),
  track_ipsig(1000, -30, 30, ""),
  initial_d0(1000, -D0RNG / 100, D0RNG / 10, "mm")
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
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

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
  int sign = (jpx*xd + jpy*yd > 0.0) ? 1 : -1;
  double ip = sign*d0;
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

  float d0 = track->trkPar[TrackParam::D0];
  hists.track_d0.fill(d0);
  float z0 = track->trkPar[TrackParam::Z0];
  hists.track_z0.fill(z0);
  {
    float d0_particle = particle->trkPar[TrackParam::D0];
    // printf("d0 %f\n", d0_particle);
    hists.particle_d0.fill(d0_particle);
    float z0_particle = particle->trkPar[TrackParam::Z0];
    hists.particle_z0.fill(z0_particle);
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
  // 	    <<  " mothers " << gen->M1 << " " << gen->M2 <<std::endl;
  hists.initial_d0.fill(std::sqrt(ix*ix + iy*iy));
}


// ____________________________________________________
// main function

int main(int argc, char *argv[])
{
  gROOT->SetBatch();
  std::string out_name("test.h5");
  if (exists(out_name) ) {
    std::cerr << out_name << " exists, exiting" << std::endl;
    return 1;
  }

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(argv[1]);

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

      bool b_label = (jet->BTag & bit::B_FLAVOR);
      bool b_tag = (jet->BTag & bit::B_TAG);

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

	// get daughter of b from top (if there is one)
	// TODO: figure out what's going on here.
	//       Look in get_parent_particle for clues...
	std::deque<int> t_seq{4, 5, 6};
	std::deque<int> at_seq{-4, -5, -6};
	GenParticle* b_daut = truth::get_parent(track, bPart, t_seq, 2);
	GenParticle* ab_daut = truth::get_parent(track, bPart, at_seq, 2);
	GenParticle* daut = b_daut ? b_daut : ab_daut;
	if (daut){
	  using std::pow;
	  using std::sqrt;
	  b_decay_pids[daut->PID]++;
	  // std::cout << truth::map_particle(daut->PID) << " "
	  // 	    << daut->PID << " "
	  // 	    << sqrt(pow(daut->X,2) + pow(daut->Y,2))
	  // 	    << std::endl;
	}

	// fill hists
	fill_track_hists(hists, track, jet);
	if (b_label) {
	  fill_track_hists(b_jet_hists, track, jet);
	} else {
	  fill_track_hists(light_jet_hists, track, jet);
	}
      }	// end loop over jet tracks
      hists.n_tracks.fill(n_tracks);
      if (track_by_pt.size() > 0) {
      	const auto* track = track_by_pt.crbegin()->second;
      	if (b_label) fill_track_hists(leading_track_b_hists, track, jet);
      	else fill_track_hists(leading_track_light_hists, track, jet);
      }
    } // end loop over jets
  }   // end loop over events
  H5::H5File out_file(out_name, H5F_ACC_EXCL);
  hists.save(out_file, "all_jets");
  b_jet_hists.save(out_file, "b_jets");
  light_jet_hists.save(out_file,"light_jets");
  leading_track_light_hists.save(out_file, "leading_track_light");
  leading_track_b_hists.save(out_file, "leading_track_b");

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
