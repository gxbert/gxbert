//
// @File: GXHadNucl3BodyMomDst.hh
//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states
//
// 20180711  Guilherme Lima -- Created, based on M.Kelsey's G4HadNucl3BodyMomDst class
//
#ifndef GXHadNucl3BodyMomDst_h
#define GXHadNucl3BodyMomDst_h 1

#include "GXVMultiBodyMomDst.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

class GXHadNucl3BodyMomDst : public GXVMultiBodyMomDst {
public:
  GXHadNucl3BodyMomDst(int verbose = 0);
  virtual ~GXHadNucl3BodyMomDst() { }
};

  /// Constructor passes arrays to templated base class
  // GXHadNucl3BodyMomDst::GXHadNucl3BodyMomDst(int verbose)
  //   : GXVMultiBodyMomDst("GXHadNucl3BodyMomDst", pqprC, psC, verbose) { }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXHadNucl3BodyMomDst_h
