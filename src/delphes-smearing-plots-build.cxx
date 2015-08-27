#include "track_set_macros.hh"
#include "root.hh"
#include "AllPlanes.hh"

#include "TFile.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include "misc_func.hh" 	// cli

#include "H5Cpp.h"
#include "ndhist/Axis.hh"
#include "ndhist/Histogram.hh"

#include <Eigen/LU>
#include <Eigen/Dense>

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
  // tracking based
  AllPlanes raveVsDelphesRaw;
  AllPlanes delphes;
};

const size_t BINS = 100;

namespace ax {
  Axis d0{"d0", BINS, -0.1, 0.1, "mm"};
  Axis z0{"z0", BINS, -1, 1, "mm"};
  Axis phi{"phi", BINS, -0.001, 0.001, "rad"};
  Axis qoverp{"qoverp", BINS, -0.005, 0.005, "GeV^{-1}"};
  Axis theta{"theta", BINS, -0.002, 0.002, "rad"};
}
typedef std::vector<Axis> Axes;
const Axes all_track_pars = {ax::d0, ax::z0, ax::phi, ax::qoverp, ax::theta};

Hists::Hists():
  raveVsDelphesRaw(all_track_pars, true),
  delphes(all_track_pars, true)
{
}

void Hists::save(H5::CommonFG& out_h5) {
#define WRITE(VAR) VAR.save_to(out_h5, #VAR)
  WRITE(raveVsDelphesRaw);
  WRITE(delphes);
#undef WRITE
}
void Hists::save(H5::CommonFG& out_file, const std::string& name) {
  H5::Group group(out_file.createGroup(name));
  save(group);
}

namespace {
  std::map<std::string, double> get_smearing(const Track& track) {
    const float* pars = track.trkPar;
    const auto* rawtrk = root::as<Track>(track.Particle.GetObject());
    const float* rawpars = rawtrk->trkPar;
    float delta[5];
    for (int iii = 0; iii < 5; iii++) {
      delta[iii] = rawpars[iii] - pars[iii];
    }
    using namespace TrackParam;
    using namespace ax;
    return {
      {d0.name, delta[D0]},
      {z0.name, delta[Z0]},
      {phi.name, delta[PHI]},
      {qoverp.name, delta[QOVERP]},
      {theta.name, delta[THETA]}
    };
  }

  typedef Eigen::Matrix<float, 5, 1> TrackParameters;
  typedef Eigen::Matrix<float, 5, 5> CovMatrix;

  TrackParameters eigenFromTrkArray(const float* trk_array) {
    TrackParameters pars;
    for (int iii = 0; iii < 5; iii++) {
      pars(iii) = trk_array[iii];
    }
    return pars;
  }

  CovMatrix eigenFromCovArray(const float* cov_array) {
    CovMatrix cov;
    // ugly copy of covariance matrix to eigen matrix
    using namespace TrackParam;
    TRKCOV_2FROMARRAY(D0);

    TRKCOV_FROMARRAY(Z0,D0);
    TRKCOV_2FROMARRAY(Z0);

    TRKCOV_FROMARRAY(PHI,D0);
    TRKCOV_FROMARRAY(PHI,Z0);
    TRKCOV_2FROMARRAY(PHI);

    TRKCOV_FROMARRAY(THETA,D0);
    TRKCOV_FROMARRAY(THETA,Z0);
    TRKCOV_FROMARRAY(THETA,PHI);
    TRKCOV_2FROMARRAY(THETA);

    TRKCOV_FROMARRAY(QOVERP,D0);
    TRKCOV_FROMARRAY(QOVERP,Z0);
    TRKCOV_FROMARRAY(QOVERP,PHI);
    TRKCOV_FROMARRAY(QOVERP,THETA);
    TRKCOV_2FROMARRAY(QOVERP);
    return cov;
  }

  double get_mahalanobis_smearing_dist(const Track& track) {
    CovMatrix cov = eigenFromCovArray(track.trkCov);
    // invert the matrix
    CovMatrix inverse = cov.inverse();
    std::cout << cov << std::endl;
    std::cout << inverse << std::endl;
    // std::cout << cov * inverse << std::endl;

    TrackParameters smeared = eigenFromTrkArray(track.trkPar);
    const auto* rawtrk = root::as<Track>(track.Particle.GetObject());
    TrackParameters raw = eigenFromTrkArray(rawtrk->trkPar);
    TrackParameters smearing = smeared - raw;
    std::cout << smearing << std::endl;
    double dist2 = smearing.transpose() * inverse * smearing;
    std::cout << dist2 << std::endl;
    double cart2 = smearing.transpose() * smearing;
    std::cout << std::sqrt(dist2) << " " << std::sqrt(cart2) << std::endl;
    abort();
    return 1;
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
  // TClonesArray* bJets = treeReader->UseBranch("Jet");
  TClonesArray* tracks_branch = treeReader->UseBranch("Track");
  treeReader->UseBranch("OriginalTrack");

  Hists hists;

  // Loop over all events
  std::cout << "looping over " << numberOfEntries << " entries" << std::endl;
  int onem = std::max(numberOfEntries / 1000, 1ll);
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    if (entry % onem == 0) std::cout << entry << " processed\r" << std::flush;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int n_tracks = tracks_branch->GetEntriesFast();

    for (int i_trk = 0; i_trk < n_tracks; i_trk++) {
      auto* track = root::as<Track>(tracks_branch->At(i_trk));
      if (std::abs(track->Eta) > 2.5) continue;
      get_mahalanobis_smearing_dist(*track);
      const auto smearing = get_smearing(*track);
      hists.delphes.fill(smearing);

    } // end loop over jets
  }   // end loop over events
  std::cout << std::endl;

  H5::H5File out_file(cli.out_name, H5F_ACC_EXCL);
  hists.save(out_file, "all");
  return 0;

}

