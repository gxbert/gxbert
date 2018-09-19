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

  static const GXMultiBodyMomentumDist* GetInstance()
  {
    if (!theInstance) {
      theInstance = new GXMultiBodyMomentumDist;
      //TK for GXBERT: Exclude to use G4AutoDelete
      //G4AutoDelete::Register(theInstance);
    }

    return theInstance;
  }

  // Return appropriate generator for initial state and multiplicity
  template <typename T>
  static const GXVMultiBodyMomDst* GetDist(Index_v<T> const& is, Index_v<T> const& mult)
  {
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

  G4ThreadLocal GXMultiBodyMomentumDist* GXMultiBodyMomentumDist::theInstance = 0;

GXMultiBodyMomentumDist::GXMultiBodyMomentumDist()
  : nn3BodyDst(new GXNuclNucl3BodyMomDst),
    nn4BodyDst(new GXNuclNucl4BodyMomDst),
    hn3BodyDst(new GXHadNucl3BodyMomDst),
    hn4BodyDst(new GXHadNucl4BodyMomDst)
{ }

void GXMultiBodyMomentumDist::passVerbose(int vb)
{
  if (nn3BodyDst) nn3BodyDst->setVerboseLevel(vb);
  if (nn4BodyDst) nn4BodyDst->setVerboseLevel(vb);
  if (hn3BodyDst) hn3BodyDst->setVerboseLevel(vb);
  if (hn4BodyDst) hn4BodyDst->setVerboseLevel(vb);
}

GXVMultiBodyMomDst const*
GXMultiBodyMomentumDist::ChooseDist(int is, int mult) const
{
  if (is == pro*pro || is == pro*neu || is == neu*neu) {
    //***** REMOVED BY VLADIMIR UZHINSKY 18 JULY 2011
    if (GXCascadeParameters::use3BodyMom() && mult==3) return nn3BodyDst;
    return nn4BodyDst;
  }

  else {	// FIXME:  All other initial states use pi-N scattering
    //***** REMOVED BY VLADIMIR UZHINSKY 18 JULY 2011
    if (GXCascadeParameters::use3BodyMom() && mult==3) return hn3BodyDst;
    return hn4BodyDst;
  }

  // Invalid interaction
  return 0;
}

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXMultiBodyMomentumDist_h
