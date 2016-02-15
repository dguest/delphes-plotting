#include "hl_var_map.hh"

#include "classes/flavortag/SecondaryVertex.hh"
#include "classes/DelphesClasses.h"

std::map<std::string, double> hl_var_map(const Jet& jet) {
  std::map<std::string, double> out;
  out["pt"] = jet.PT;
  out["eta"] = jet.Eta;
  out["track_2_d0_significance"] = jet.track2d0sig;
  out["track_3_d0_significance"] = jet.track3d0sig;
  out["track_2_z0_significance"] = jet.track2z0sig;
  out["track_3_z0_significance"] = jet.track3z0sig;
  out["n_tracks_over_d0_threshold"] = jet.tracksOverIpThreshold;
  out["jet_prob"] = jet.jetProb;
  out["jet_width_eta"] = jet.jetWidthEta;
  out["jet_width_phi"] = jet.jetWidthPhi;
  const THighLevelSecondaryVertex& svx = jet.HLSecondaryVertex;
  out["vertex_significance"] = svx.svLsig;
  out["n_secondary_vertices"] = svx.svNVertex;
  out["n_secondary_vertex_tracks"] = svx.svNTracks;
  out["delta_r_vertex"] = svx.svDrJet;
  out["vertex_mass"] = svx.svMass;
  out["vertex_energy_fraction"] = svx.svEnergyFraction;
  return out;
}

