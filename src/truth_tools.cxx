#include "truth_tools.hh"

#include "TClonesArray.h"
#include "TClass.h"

#include "classes/DelphesClasses.h"
#include "root.hh"

#include <deque>

// this stuff should normally be silent
#ifdef DEBUG
#include <iostream>
#endif // DEBUG

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

  int walk_idx(TClonesArray* particles, int idx, int target = 25) {
    if (idx == -1) return -1;
    GenParticle* part = root::as<GenParticle>(particles->At(idx));
    // printf("PID: %i\n", part->PID);
    if (std::abs(part->PID) == target) return idx;
    int first_try = walk_idx(particles, part->M1, target);
    if (first_try != -1) return first_try;
    return walk_idx(particles, part->M2, target);
  }

  typedef std::deque<int> ISEQ;
  ISEQ walk_pids(TClonesArray* particles, int idx,
		 const ISEQ& target = {5, 6}, ISEQ history = {0, 0},
		 ISEQ indices = {-2, -2}) {
    indices.push_back(idx);
    indices.pop_front();
    if (idx == -1) return indices;
    GenParticle* part = root::as<GenParticle>(particles->At(idx));
    int pid = part->PID;
    history.push_back(std::abs(pid));
    history.pop_front();

    if (history == target){
      return indices;
    }

    ISEQ first_try = walk_pids(particles, part->M1, target, history, indices);
    if (first_try.back() != -1) {
      // printf(" %i", history.back());
      first_try.push_front(idx);
      return first_try;
    }
    ISEQ second_try = walk_pids(particles, part->M2, target, history, indices);
    if (second_try.back() != -1) {
      // printf(" %i", history.back());
      second_try.push_front(idx);
    }
    return second_try;
  }

  GenParticle* get_daughter(TClonesArray* particles, int idx,
			    const ISEQ target, bool verb = false) {
    ISEQ seq = walk_pids(particles, idx, target);
    int mom_idx = seq.back();
    if (mom_idx == -1) return 0;
    int daughter_idx = seq.at(seq.size() - target.size());

#ifdef DEBUG
    // IODEBUG
    std::cout << "found!";
    for (auto idx: seq) {
      std::cout << " " << root::as<GenParticle>(particles->At(idx))->PID;
    }
    std::cout << std::endl;
#endif //DEBUG

    return root::as<GenParticle>(particles->At(daughter_idx));
  }

}

namespace truth {

  GenParticle* get_gen_particle(const Track* track) {
    int depth = 0;
    TObject* thing = walk_track(track, depth);
    if (thing) return root::as<GenParticle>(thing);
    return 0;
  }



  GenParticle* get_parent_with_decay(const Track* track,
				     TClonesArray* particles,
				     ISEQ sequence){
    GenParticle* part = get_gen_particle(track);
    if (!part) return 0;
    GenParticle* decay1 = get_daughter(particles, part->M1, sequence);
    if (decay1) return decay1;
    return get_daughter(particles, part->M2, sequence);
  }

}
