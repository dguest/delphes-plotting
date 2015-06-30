#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TSystem.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"

#include "ExRootTreeReader.h"
#include "misc_func.hh"

#include "H5Cpp.h"
#include "Histogram.hh"

const unsigned MAX_TRACKS = 200;
const double GeV = 1;
const double MeV = 0.001;

struct Hists {
  Hists();
  void save(std::string);
  Histogram n_tracks;
  Histogram track_pt;
  Histogram track_d0;
  Histogram particle_d0;
  Histogram track_z0;
  Histogram particle_z0;
  Histogram track_d0sig;
  Histogram track_z0sig;
};

const float D0RNG = 0.1;
const float Z0RNG = 0.1;

Hists::Hists():
  n_tracks(MAX_TRACKS, -0.5, MAX_TRACKS + 0.5),
  track_pt(200, 0, 200, "GeV"),
  track_d0(100, -D0RNG, D0RNG, "mm"),
  particle_d0(100, -D0RNG, D0RNG, "mm"),
  track_z0(100, -Z0RNG, Z0RNG, "mm"),
  particle_z0(100, -Z0RNG, Z0RNG, "mm"),
  track_d0sig(100, -10, 10),
  track_z0sig(100, -10, 10)
{
}
void Hists::save(std::string output) {
  H5::H5File out_file(output, H5F_ACC_EXCL);
#define WRITE(VAR) VAR.write_to(out_file, #VAR)
  WRITE(n_tracks);
  WRITE(track_pt);
  WRITE(track_d0);
  WRITE(particle_d0);
  WRITE(track_z0);
  WRITE(particle_z0);
  WRITE(track_d0sig);
  WRITE(track_z0sig);
#undef WRITE
}


int main(int argc, char *argv[])
{
  gROOT->SetBatch();
  // absolute path to delphes has to be passed in as a macro DELPHESDIR
  gSystem->Load(DELPHESDIR "/libDelphes");
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
  TClonesArray* bTrack = treeReader->UseBranch("Track");
  TClonesArray* bOriginalTrack = treeReader->UseBranch("OriginalTrack");

  Hists hists;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int n_tracks = bTrack->GetEntries();
    hists.n_tracks.fill(n_tracks);
    for (int i_track = 0; i_track < n_tracks; i_track++) {
      Track* track = (Track*) bTrack->At(i_track);

      if (track->PT < 1*GeV) continue;

      hists.track_pt.fill(track->PT);
      float d0 = track->trkPar[TrackParam::D0];
      hists.track_d0.fill(d0);
      float z0 = track->trkPar[TrackParam::Z0];
      hists.track_z0.fill(z0);
      {
	Track* particle = (Track*) track->Particle.GetObject();
	float d0_particle = particle->trkPar[TrackParam::D0];
	hists.particle_d0.fill(d0_particle);
	float z0_particle = particle->trkPar[TrackParam::Z0];
	hists.particle_z0.fill(z0_particle);
      }
      {
	float d0sig = d0 / std::sqrt(track->trkCov[TrackParam::D0D0]);
	hists.track_d0sig.fill(d0sig);
      }
      {
	float z0sig = z0 / std::sqrt(track->trkCov[TrackParam::Z0Z0]);
	hists.track_z0sig.fill(z0sig);
      }
    }

  }
  hists.save(out_name);
  return 0;

}
