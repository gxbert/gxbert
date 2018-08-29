//
// @File: GXNuclNucl4BodyMomDst.hh
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states
//
// 20180711  Guilherme Lima -- Created, based on M.Kelsey's G4NuclNucl4BodyMomDst class
//
#ifndef GXNuclNucl4BodyMomDst_h
#define GXNuclNucl4BodyMomDst_h 1

#include "GXVMultiBodyMomDst.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

class GXNuclNucl4BodyMomDst : public GXVMultiBodyMomDst {
public:
  GXNuclNucl4BodyMomDst(int verbose = 0);
  virtual ~GXNuclNucl4BodyMomDst() { }
};

  /// Constructor passes arrays to templated base class
  // GXNuclNucl4BodyMomDst::GXNuclNucl4BodyMomDst(int verbose)
  //   : GXVMultiBodyMomDst("GXNuclNucl4BodyMomDst", pqprC, psC, verbose) { }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXNuclNucl4BodyMomDst_h
