#include "truth_tools.hh"

#include "TClonesArray.h"
#include "TClass.h"

#include "classes/DelphesClasses.h"
#include "root.hh"

#include <deque>
#include <cassert>
#include <iostream>
#include <set>

namespace {

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
  ISEQ walk_pids(TClonesArray* particles, int idx,
		 const ISEQ& target = {5, 6}, ISEQ history = {0, 0}) {

    // return error code if we were directed to an empty particle
    if (idx == -1){
      return {-1};
    }

    // get the particle, check if it satisfies the sequence
    GenParticle* part = root::as<GenParticle>(particles->At(idx));
    int pid = std::abs(part->PID); // NOTE: do this to ignore anti-particles
    // we only pay attention when the PID changes
    if (pid != history.back()){
      history.push_back(pid);
    }
    if (history.size() > target.size()) history.pop_front();

    if (history == target){
      return {idx};
    }

    for (int nextidx: {part->M1, part->M2} ) {
      ISEQ seq = walk_pids(particles, nextidx, target, history);
      if (seq.back() != -1) {
	seq.push_front(idx);
	return seq;
      }
    }
    return {-1};
  }

  // some PIDs are ignored while walking the decay chain
  std::set<int> ignored_pid({91, 92, 93});

  GenParticle* get_parent_particle(TClonesArray* particles, int start_idx,
				   const ISEQ target, int step_back) {
    using root::as;
    ISEQ seq = walk_pids(particles, start_idx, target);
    int mom_idx = seq.back();
    if (mom_idx == -1) return 0;

    int last_pid = std::abs(as<GenParticle>(particles->At(mom_idx))->PID);
    assert(last_pid == target.back());

    // walk back some number of steps to find the daughter of the decay
    for (auto idx = seq.rbegin(); idx != seq.rend(); idx++) {
      GenParticle* part = as<GenParticle>(particles->At(*idx));
      if (step_back == 0){
	return part;
      }
      int this_pid = std::abs(part->PID);
      if (!ignored_pid.count(this_pid) && last_pid != this_pid) {
	last_pid = this_pid;
	step_back--;
      }
    }
    return 0;
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
      GenParticle* decay = get_parent_particle(
	particles, id, sequence, step_back);
      if (decay) return decay;
    }
    return 0;
  }

}
