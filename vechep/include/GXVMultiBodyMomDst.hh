//
// @File: GXVMultiBodyMomDst.hh
//
// Description: intermediate base class for INUCL parametrizations of
//		final-state momentum distributions in Bertini-style cascade
//
// NOTE:  Coefficient arrays have fixed dimensions, validated by compiler
//
// 20180711  Guilherme Lima -- Created, combining M.Kelsey's G4VMultiBodyMomDst and G4InuclParamMomDst classes
//
#ifndef GXVMultiBodyMomDst_h
#define GXVMultiBodyMomDst_h 1

#include "VecCore/VecCore"
#include "G4InuclParticleNames.hh"
#include "GXPowVec.hh"
#include <string>

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

using namespace G4InuclParticleNames;
  
class GXVMultiBodyMomDst {
public:
  // NOTE:  Array arguments must be STATIC, GLOBAL declarations
  GXVMultiBodyMomDst(const std::string& name, 
		     const double (&pqprC)[2][4][4],
		     const double (&psC)[2][3],
		     int verbose=0)
    : fVerbose(verbose)
    , fName(name)
    , coeffPR(pqprC)
    , coeffPS(psC)
  { }

  virtual ~GXVMultiBodyMomDst() { }
  
  template <typename T, typename IntT>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T GetMomentum(IntT const& ptype, T const& ekin) const;

  void setVerboseLevel(int vb=0) { fVerbose = vb; }
  const std::string& GetName() const { return fName; }

protected:
  int fVerbose;
  std::string fName;
  const double (&coeffPR)[2][4][4];	// (coeffs Ekin^0..3) * S^0..3
  const double (&coeffPS)[2][3];	// PS = sum coeffs * Ekin^0..2 
};        

  // Use coefficients in power expansion of random fraction

  template <typename T, typename IntT>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T GXVMultiBodyMomDst::GetMomentum(IntT const& ptype, T const& ekin) const
  {
    if (fVerbose>3) {
      std::cerr << fName << "::GetMomentum: ptype=" << ptype << " ekin=" << ekin <<"\n";
    }
    if (!isHomogeneous(ptype)) {
      std::cerr <<"GXVMultiBodyMomDst::GetMomentum(ptype,ekin) called with non-homogeneous ptype: "<< ptype <<"\n";
      return T(0.0);
    }

    int itype = Get(ptype,0);
    int JK = (itype==pro || itype==neu) ? 0 : 1;	// nucleon vs. other

    if (fVerbose > 3) std::cerr << " JK " << JK <<"\n";

    GXPowVec<T,IntT>* theGXPow = GXPowVec<T,IntT>::GetInstance();	// For convenience
    T Spow = randomInuclPowers(ekin, coeffPR[JK]);

    double C=0.;
    T PS(0.);
    for(int im = 0; im < 3; im++) {
      C = coeffPS[JK][im];
      PS += C * theGXPow->PowN(ekin, im);

      if (fVerbose >3) {
	std::cerr << " im " << im << " : coeffPS[JK][im] " << C
	     << " ekin^im " << theGXPow->PowN(ekin, im) <<"\n";
      }
    }
  
    T PRA = PS * Spow;
    if (fVerbose > 3) {
      std::cerr << " PS=" << PS << " Spow = sqrt(S)*(PR+(1-PQ)*S^4) = " << Spow
	     << " and PRA = PS*Spow = " << PRA <<"\n";
    }

    return math::Abs(PRA);
  }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXVMultiBodyMomDst_h
