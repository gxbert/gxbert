//
// @File: GXHadNucl4BodyMomDst.hh
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states
//
// 20180711  Guilherme Lima -- Created, based on M.Kelsey's G4HadNucl4BodyMomDst class
//
#ifndef GXHadNucl4BodyMomDst_h
#define GXHadNucl4BodyMomDst_h 1

#include "GXVMultiBodyMomDst.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

class GXHadNucl4BodyMomDst : public GXVMultiBodyMomDst {
public:
  GXHadNucl4BodyMomDst(int verbose = 0);
  virtual ~GXHadNucl4BodyMomDst() { }
};

  /// Constructor passes arrays to templated base class
  // GXHadNucl4BodyMomDst::GXHadNucl4BodyMomDst(int verbose)
  //   : GXVMultiBodyMomDst("GXHadNucl4BodyMomDst", pqprC, psC, verbose) { }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXHadNucl4BodyMomDst_h
