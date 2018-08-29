//
// @File: GXMultiBodyMomentumDist.hh
//
// Description: Singleton class to evaluate multi-body momentum distribution
//		functions based on initial state codes and multiplicity.
//
// 2018-08-15  Guilherme Lima  - Created, based on M.Kelsey's G4MultiBodyMomentumDist class

#ifndef GXMultiBodyMomentumDist_h
#define GXMultiBodyMomentumDist_h 1

#include "GXNuclNucl3BodyMomDst.hh"
#include "GXNuclNucl4BodyMomDst.hh"
#include "GXHadNucl3BodyMomDst.hh"
#include "GXHadNucl4BodyMomDst.hh"

class GXVMultiBodyMomDst;

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {


class GXMultiBodyMomentumDist {
public:
  ~GXMultiBodyMomentumDist();

  static const GXMultiBodyMomentumDist* GetInstance();

  // Return appropriate generator for initial state and multiplicity
  static const GXVMultiBodyMomDst* GetDist(Int_v is, Int_v mult) {
    if( !isHomogeneous(is) ) return NULL;
    if( !isHomogeneous(mult) ) return NULL;
    return GetInstance()->ChooseDist(Get(is,0), Get(mult,0));
  }

  // Pass verbosity through to owned objects
  static void setVerboseLevel(int vb=0) {
    const_cast<GXMultiBodyMomentumDist*>(GetInstance())->passVerbose(vb);
  }

private:
  // Constructor is private for singleton
  GXMultiBodyMomentumDist();

  GXVMultiBodyMomDst const* ChooseDist(int is, int mult) const;

  void passVerbose(int verbose);	// Pass verbosity through instance

  static G4ThreadLocal GXMultiBodyMomentumDist* theInstance;	// Per thread

  // Generators for various initial/final state combinations
  GXNuclNucl3BodyMomDst* nn3BodyDst;	// N N to X Y Z
  GXNuclNucl4BodyMomDst* nn4BodyDst;    // N N to X Y Z W ...
  GXHadNucl3BodyMomDst*  hn3BodyDst;	// pi,K,Y,gamma N to X Y Z
  GXHadNucl4BodyMomDst*  hn4BodyDst;	// pi,K,Y,gamma N to X Y Z W ...

private:
  // Copying of modules is forbidden
  GXMultiBodyMomentumDist(const GXMultiBodyMomentumDist&);
  GXMultiBodyMomentumDist& operator=(const GXMultiBodyMomentumDist&);
};

  void GXMultiBodyMomentumDist::passVerbose(int vb)
  {
    if (nn3BodyDst) nn3BodyDst->setVerboseLevel(vb);
    if (nn4BodyDst) nn4BodyDst->setVerboseLevel(vb);
    if (hn3BodyDst) hn3BodyDst->setVerboseLevel(vb);
    if (hn4BodyDst) hn4BodyDst->setVerboseLevel(vb);
  }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXMultiBodyMomentumDist_h
