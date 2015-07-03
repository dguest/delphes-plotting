#ifndef TRUTH_TOOLS
#define TRUTH_TOOLS

class GenParticle;
class Track;
class TClonesArray;

#include <deque>

namespace truth {
  GenParticle* get_parent_with_decay(
    const Track* track, TClonesArray* particles,
    std::deque<int> sequence = {5, 6});

  GenParticle* get_gen_particle(const Track* track);
}

#endif
