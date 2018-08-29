//
// @File: GXNuclNucl3BodyMomDst.hh
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states
//
// 20180711  Guilherme Lima -- Created, based on M.Kelsey's G4NuclNucl3BodyMomDst class
//
#ifndef GXNuclNucl3BodyMomDst_h
#define GXNuclNucl3BodyMomDst_h 1

#include "GXVMultiBodyMomDst.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

class GXNuclNucl3BodyMomDst : public GXVMultiBodyMomDst {
public:
  GXNuclNucl3BodyMomDst(int verbose = 0);
  virtual ~GXNuclNucl3BodyMomDst() { }
};

  /// Constructor passes arrays to templated base class
  // GXNuclNucl3BodyMomDst::GXNuclNucl3BodyMomDst(int verbose)
  //   : GXVMultiBodyMomDst("GXNuclNucl3BodyMomDst", pqprC, psC, verbose) { }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXNuclNucl3BodyMomDst_h
