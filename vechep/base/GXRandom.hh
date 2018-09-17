#ifndef GXRANDOM_H
#define GXRANDOM_H

#include "VecCoreLib/Rng/MRG32k3a.h"
#include "VecHepDefs.hh"

namespace gxbert {
inline namespace GXBERT_IMPL_NAMESPACE {

class GXRandom 
{
public:
  static vecRng::MRG32k3a<VectorBackend>& RNG();

public:

  VECCORE_ATT_HOST 
  GXRandom() {}

  VECCORE_ATT_HOST
  ~GXRandom() {}

  //add utilities
  //redesign for multithreaded applications
};

vecRng::MRG32k3a<VectorBackend>& GXRandom::RNG(){
  static vecRng::cxx::MRG32k3a<VectorBackend> theGenerator 
    = vecRng::cxx::MRG32k3a<VectorBackend>();
  theGenerator.Initialize();
  return theGenerator;
}

} // end namespace impl
} // end namespace gxbert

#endif
