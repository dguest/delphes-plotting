#include <iostream>
#include <utility>
#include <vector>
#include <string>

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

struct Hists {
  Hists();
  void save(std::string);
  Histogram n_tracks;
  Histogram track_d0;
  Histogram particle_d0;
};
Hists::Hists():
  n_tracks(MAX_TRACKS, -0.5, MAX_TRACKS + 0.5),
  track_d0(100, -0.25, 0.25),
  particle_d0(100, -0.25, 0.25)
{
}
void Hists::save(std::string output) {
  H5::H5File out_file(output, H5F_ACC_EXCL);
  n_tracks.write_to(out_file, "n_tracks");
  track_d0.write_to(out_file, "track_d0");
  particle_d0.write_to(out_file, "particle_d0");
}

int main(int argc, char *argv[])
{
  gROOT->SetBatch();
  gSystem->Load("delphes/libDelphes");
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
      float d0 = track->trkPar[TrackParam::D0];
      hists.track_d0.fill(d0);
      Track* particle = (Track*) track->Particle.GetObject();
      float d0_particle = particle->trkPar[TrackParam::D0];
      hists.particle_d0.fill(d0_particle);
    }

  }
  hists.save(out_name);
  return 0;

}
