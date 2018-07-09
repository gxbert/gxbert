//
// @File: GXParticleLargerEkin.h
//
// 20180627  Guilherme Lima -- Created, based on M.Kelsey's G4ParticleLargerEkin class
//
#ifndef GXPARTICLE_LARGEREKIN_H
#define GXPARTICLE_LARGEREKIN_H

//#include "G4CascadParticle.hh"
#include "GXInuclElementaryParticle.hh"

#ifdef G4CASCADE_DEBUG_SORT
  #include <iostream>
#endif

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

template <typename T>
class GXParticleLargerEkin {

  using Bool_v = vecCore::Mask_v<T>;

public:
  Bool_v operator() (const GXInuclElementaryParticle<T>& part1,
		     const GXInuclElementaryParticle<T>& part2) {
#ifdef G4CASCADE_DEBUG_SORT
    cerr << "part1 @ " << &part1 << ": " << part1
	 << "part2 @ " << &part2 << ": " << part2 <<"\n";
#endif
    return (part1.getKineticEnergy() > part2.getKineticEnergy());
  }
 
  Bool_v operator() (const GXInuclElementaryParticle<T>* part1,
		     const GXInuclElementaryParticle<T>* part2) {
    return (part1 && part2 && operator()(*part1, *part2));
  }

  /*
  Bool_V operator() (const G4CascadParticle& part1,
		     const G4CascadParticle& part2) {
    return (operator()(part1.getParticle(), part2.getParticle()));
  }

  Bool_V operator() (const G4CascadParticle* part1,
		     const G4CascadParticle* part2) {
    return (part1 && part2 && operator()(*part1, *part2));
  }
  */
};

}
}
#endif // GXPARTICLE_LARGEREKIN_H
