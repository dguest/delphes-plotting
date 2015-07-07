#include "truth_tools.hh"

#include "TClonesArray.h"
#include "TClass.h"

#include "classes/DelphesClasses.h"
#include "root.hh"

#include <deque>
#include <cassert>
#include <iostream>
#include <set>
#include <map>
#include <list>

namespace {
  // we ignore `minor` particle changes
  bool is_significant_shift(int pid1, int pid2);
  int major_quark(int pid);
  bool is_metastable_hadron_parent(int pid);

  TObject* walk_track(const Track* track, int& depth){
    if (!track) return 0;
    TObject* particle = track->Particle.GetObject();
    if (!particle) return 0;
    if (root::is<Track>(particle)){
      depth++;
      return walk_track(root::as<Track>(particle), depth);
    }
    if (root::is<GenParticle>(particle)) {
      return particle;
    }
    printf("fuck!\n");
    return 0;
  }

  /* not currently used for anything
  int walk_idx(TClonesArray* particles, int idx, int target = 25) {
    if (idx == -1) return -1;
    GenParticle* part = root::as<GenParticle>(particles->At(idx));
    // printf("PID: %i\n", part->PID);
    if (std::abs(part->PID) == target) return idx;
    int first_try = walk_idx(particles, part->M1, target);
    if (first_try != -1) return first_try;
    return walk_idx(particles, part->M2, target);
  }
  */

  // main recursive truth record walk function
  typedef std::deque<int> ISEQ;
  ISEQ walk_particle_index(TClonesArray* particles, int idx,
			   const ISEQ& target = {5, 6}, ISEQ history = {-1}) {

    // return error code if we were directed to an empty particle
    if (idx == -1){
      return {-1};
    }

    // get the particle, check if it satisfies the sequence
    GenParticle* part = root::as<GenParticle>(particles->At(idx));
    int pid = part->PID;
    // we only pay attention when the PID changes
    if (is_significant_shift(pid,history.back())){
      history.push_back(std::copysign(major_quark(pid), pid));
    }
    if (history.size() > target.size()) history.pop_front();

    if (history == target){
      return {idx};
    }

    for (int nextidx: {part->M1, part->M2} ) {
      ISEQ seq = walk_particle_index(particles, nextidx, target, history);
      if (seq.back() != -1) {
	seq.push_front(idx);
	return seq;
      }
    }
    return {-1};
  }

  //  I've looked up some PIDs
  std::map<int, std::string> pid_map ({
      // most common from t decays
      {511, "B0"},
      {521, "B+"},
      {531, "Bs0"},
      {413, "D*+"},
      {5122, "Lb0"},
      {423, "D*0"},
      {311, "K0"},
      {113, "rho0"},
      {513, "B*0"},
      {523, "B*+"},
      {310, "Ks0"},
      {433, "D*s+"},
      {213, "rho+"},
      {4122, "Lc+"},
      {421, "D0"},
      {223, "omega"},
      {411, "D+"},
      {3122, "L0"},
      {431, "Ds+"},
      {221, "eta"},
      {313, "K*0"},
      {1, "d"},
      {2, "u"},
      {323, "K*+"},
      {3, "s"},
      {331, "eta'"},
      {15, "tau-"},
      {20213, "a1+"},
      {443, "J/psi"},
      {111, "pi0"},
      {415, "D*2+"},
      {425, "D*20"},
      {4132, "Xic0"},
      {441, "etac"},
	// other useful SM ones
      {4, "c"},
      {5, "b"},
      {6, "t"},
      {21, "g"},
      {22, "gamma"},
      {23, "Z"},
      {24, "W+"},
      {25, "h"}
    });
  // change the sign of anti-particles
  std::string flop_sign(std::string base) {
      auto pos = base.find("+");
      if (pos != std::string::npos) {
	base.replace(pos, 1, "-");
	return base;
      }
      pos = base.find("-");
      if (pos != std::string::npos) {
	base.replace(pos, 1, "+");
	return base;
      }
      return "anti-" + base;
  }

  // some PIDs are ignored while walking the decay chain
  std::set<int> ignored_pid({91, 92, 93});

  // get heaviest quark in the pdgid
  int major_quark(int pid){
    int absid = std::abs(pid);
    if (absid <= 6) return absid;
    int tens = (absid % 100) / 10;
    int hundreds = (absid % 1000) / 100;
    int thousands = (absid % 10000) / 1000;
    // things with no quarks (tens and thousands == 0) just get abs pid
    if (hundreds == 0 && thousands == 0) {
      return absid;
    }
    // the leading digit should be larger than or equal to the lower ones
    // special exception for K_short, number 130
    assert(tens <= hundreds || tens <= thousands || absid == 130);
    assert(thousands == 0 || hundreds <= thousands);
    return std::max({tens, hundreds, thousands});
  }

  bool is_significant_shift(int pid1, int pid2) {
    int absid1 = std::abs(pid1);
    int absid2 = std::abs(pid2);

    // index with `major_quark` function above
    int major1 = major_quark(absid1);
    int major2 = major_quark(absid2);
    return (major1 != major2);
  }

  // get hadron number as q1q2q3
  int hadron_number(int pid) {
    int absid = std::abs(pid);
    if (absid < 100) return 0;
    return (absid / 10) % 1000;
  }

  // check for particle that decays to a weekly-decaying hadron
  bool is_metastable_hadron_parent(int pid) {
    int absid = std::abs(pid);
    // all bare quarks below t should hadronize
    if (absid <= 5) return true;

    // hadron number is of the form q1q2q3
    int had_number = hadron_number(absid);
    if (had_number == 0) return false;

    // things with high spin decay
    int spin = absid % 10;
    if (spin > 2) return true;
    // anything left is a metastable hadron
    return false;
  }

  // remove `insignificant' particles from the decay chain. When two
  // particles differ insignificantly, the upstream one is removed.
  // Originally implemented as a list (thus the complexity).
  void clean_particles(std::vector<GenParticle*>& parts){
    if (parts.size() < 2) return;
    auto iter1 = parts.begin();
    auto iter2 = iter1;
    iter2++;
    while(iter2 != parts.end()) {
      const auto& part1 = **iter1;
      const auto& part2 = **iter2;
      bool ignored = ignored_pid.count(part2.PID);
      bool insig = !is_significant_shift(part1.PID, part2.PID);
      bool had_parent = is_metastable_hadron_parent(part2.PID);
      if (ignored || insig || had_parent) {
	iter2 = parts.erase(iter2);
      } else {
	iter2++;
      }
      iter1 = iter2;
      iter1--;
    }
  }

  GenParticle* get_chain(TClonesArray* particles, int start_idx,
			 const ISEQ target, int step_back) {
    using root::as;
    std::deque<int> end_target{target.back()};
    ISEQ seq = walk_particle_index(particles, start_idx, end_target);
    int mom_idx = seq.back();
    if (mom_idx == -1) return 0;

    assert(seq.front() == start_idx);

    // grab every particle along the path
    std::vector<GenParticle*> parts;
    for (const auto& idx: seq) {
      parts.push_back(as<GenParticle>(particles->At(idx)));
    }
    // remove `insignificant' decays
    clean_particles(parts);

    // check to see if particle we're asking for is out of range...
    if (target.size() > parts.size() || step_back >= parts.size() ) {
      return 0;
    }
    int pend = parts.size() - 1;
    int tend = target.size() - 1;
    for (int iii = 0; iii <= tend; iii++) {
      int targ = target.at(tend - iii);
      const auto& part = *parts.at(pend - iii);
      bool sig_shift = is_significant_shift(targ, part.PID);
      bool wrong_sign = std::signbit(targ) != std::signbit(part.PID);
      // bool wrong_sign = false;
      if (sig_shift || wrong_sign) return 0;
    }

    return parts.at(pend - step_back);
  }



}

namespace truth {

  GenParticle* get_gen_particle(const Track* track) {
    int depth = 0;
    TObject* thing = walk_track(track, depth);
    if (thing) return root::as<GenParticle>(thing);
    return 0;
  }


  GenParticle* get_parent(const Track* track,
			  TClonesArray* particles,
			  ISEQ sequence, int step_back){
    GenParticle* part = get_gen_particle(track);
    if (!part) return 0;
    for (int id: {part->M1, part->M2} ){
      GenParticle* decay = get_chain(
	particles, id, sequence, step_back);
      if (decay) return decay;
    }
    return 0;
  }

  std::string map_particle(int pid) {
    bool neg = std::signbit(pid);
    int apid = std::abs(pid);
    if (!pid_map.count(apid)) return std::to_string(pid);
    auto base = pid_map.at(std::abs(pid));
    if (neg) base = flop_sign(base);
    return base;
  }
}
