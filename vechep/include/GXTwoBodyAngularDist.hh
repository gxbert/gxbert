//
// @File:    GXTwoBodyAngularist.hh
//
// Description: Singleton class to evaluate two-body angular distribution
//		functions based on intial/final state codes.
//
// 2018-08-15 Guilherme Lima  - Created, based on M.Kelsey's G4TwoBodyAngularDist.hh
//
#ifndef GXTwoBodyAngularDist_h
#define GXTwoBodyAngularDist_h 1

#include "G4TwoBodyAngularDist.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {
  /*
class G4VTwoBodyAngDst;
class G4GamP2NPipAngDst;
class G4GamP2PPi0AngDst;
class G4PP2PPAngDst;
class G4NP2NPAngDst;
class G4NuclNuclAngDst;
class G4Pi0P2Pi0PAngDst;
class G4PimP2Pi0NAngDst;
class G4PimP2PimPAngDst;
class G4PipP2PipPAngDst;
class G4PiNInelasticAngDst;
class G4HadNElastic1AngDst;
class G4HadNElastic2AngDst;
class G4GammaNuclAngDst;

class G4HadNucl3BodyAngDst;	// TEMPORARY, until migration to GENBOD
class G4NuclNucl3BodyAngDst;
  */

template <typename T>
class GXTwoBodyAngularDist {
public:
  constexpr static size_t fvsize = vecCore::VectorSize<T>();
  using Int_v = typename vecCore::backend::VcSimdArray<fvsize>::Int_v;

  ~GXTwoBodyAngularDist();

  static const G4TwoBodyAngularDist* GetInstance() { return G4TwoBodyAngularDist::GetInstance(); }

  // Return appropriate generator for initial, final state, and kw flag
  static void GetDist(Index_v<T> const& is, Index_v<T> const& fs, int kw, const G4VTwoBodyAngDst** angDist) {
    std::cerr<<"GXTwoBodyAngularDist::GetDist(): is="<< is <<"\n";
    assert( isHomogeneous(is) && "Non-homogeneous initial state!");

    int isi = Get(is, 0);
    for(size_t i = 0; i < fvsize; ++i) { angDist[i] = GetInstance()->GetDist(isi, Get(fs,i), kw); }
  }

  static void GetDist(Index_v<T> const& is, const G4VTwoBodyAngDst** angDist) {
    assert( !isHomogeneous(is) );
    int isi = Get(is, 0);
    for(size_t i = 0; i < fvsize; ++i) { angDist[i] = GetInstance()->GetDist(isi, 0, 0); }
  }

  static G4VTwoBodyAngDst const* GetDist(int is)
  {
    return G4TwoBodyAngularDist::GetInstance()->GetDist(is, 0, 0);
  }

  // Pass verbosity through to owned objects
  static void setVerboseLevel(int vb=0)
  {
    GetInstance()->setVerboseLevel(vb);
  }

private:
  // Constructor is private for singleton
  GXTwoBodyAngularDist();

  // const G4VTwoBodyAngDst* ChooseDist(int is, int fs, int kw) const
  // {
  //   return GetDist(is, fs, kw );
  // }

  static G4ThreadLocal GXTwoBodyAngularDist* theInstance;	// Per thread
  static G4ThreadLocal G4VTwoBodyAngDst* fAngDstInstance;

  // Generators for various initial/final state combinations
  G4GamP2NPipAngDst* gp_npip;		// gamma p -> n pi+
  G4GamP2PPi0AngDst* gp_ppi0;		// gamma p -> p pi0
  G4PP2PPAngDst* ppAngDst;              // pp, nn elastic
  G4NP2NPAngDst* npAngDst;              // np and pn elastic
  G4NuclNuclAngDst* nnAngDst;		// Y N elastic and inelastic
  G4Pi0P2Pi0PAngDst* pi0pAngDst;        // pi0 p, pi0 n elastic
  G4PimP2Pi0NAngDst* pipCXAngDst;       // pi- p, pi+ n, pi0 p, pi0 n charge exchange
  G4PimP2PimPAngDst* pimpAngDst;        // pi- p, pi+ n elastic
  G4PipP2PipPAngDst* pippAngDst;        // pi+ p, pi- n elastic

  G4PiNInelasticAngDst* qxAngDst;	// pi N charge/strangeness exchange
  G4HadNElastic1AngDst* hn1AngDst;	// pi+p and related elastic scattering
  G4HadNElastic2AngDst* hn2AngDst;	// pi-p and related elastic scattering
  G4GammaNuclAngDst* gnAngDst;		// gamma N inelastic

  // TEMPORARY generators for three-body final states
  G4HadNucl3BodyAngDst* hn3BodyDst;	// (pi,K,Y,g) N -> XYZ scattering
  G4NuclNucl3BodyAngDst* nn3BodyDst;	// N N -> XYZ scattering

private:
  // Copying of modules is forbidden
  GXTwoBodyAngularDist(const GXTwoBodyAngularDist&);
  GXTwoBodyAngularDist& operator=(const GXTwoBodyAngularDist&);
};

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif	// GXTwoBodyAngularDist_h
