//
// @File: GXInuclParamMomDst.hh
//
// Description: intermediate base class for INUCL parametrizations of
//		final-state momentum distributions in Bertini-style cascade
//
// NOTE:  Coefficient arrays have fixed dimensions, validated by compiler
//
// 20180711  Guilherme Lima -- Created, based on M.Kelsey's G4InuclParamMomDst class
//
#ifndef GXInuclParamMomDst_h
#define GXInuclParamMomDst_h 1

#include "G4VMultiBodyMomDst.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

class GXVMultiBodyMomDst {
public:
  // NOTE:  Array arguments must be STATIC, GLOBAL declarations
  G4InuclParamMomDst(const G4String& name, 
		     const G4double (&pqprC)[2][4][4],
		     const G4double (&psC)[2][3],
		     G4int verbose=0)
    : fName(name)
    , coeffPR(pqprC)
    , coeffPS(psC)
    , fVerbose(verbose)
  { }

  virtual ~G4InuclParamMomDst() { }
  
  template <typename T>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T GetMomentum(vecCore::Index_v<T> const& ptype, T const& ekin) const;

protected:
  int fVerbose;
  std::string fName;
  const G4double (&coeffPR)[2][4][4];	// (coeffs Ekin^0..3) * S^0..3
  const G4double (&coeffPS)[2][3];	// PS = sum coeffs * Ekin^0..2 
};        

  // Use coefficients in power expansion of random fraction

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T G4InuclParamMomDst::GetMomentum(G4int ptype, const G4double& ekin) const
  {
    if (verboseLevel>3) {
      cerr << theName << "::GetMomentum: ptype " << ptype << " ekin " << ekin <<"\n";
    }

    G4int JK = (ptype==pro || ptype==neu) ? 0 : 1;	// nucleon vs. other

    if (verboseLevel > 3) G4cout << " JK " << JK << G4endl;

    GXPow* theGXPow = GXPow::GetInstance();	// For convenience

    G4double Spow = randomInuclPowers(ekin, coeffPR[JK]);

    G4double C=0., PS=0.;
    for(G4int im = 0; im < 3; im++) {
      C = coeffPS[JK][im];
      PS += C * theGXPow->powN(ekin, im);

      if (verboseLevel >3) {
	G4cout << " im " << im << " : coeffPS[JK][im] " << C
	       << " ekin^im " << theGXPow->powN(ekin, im) << G4endl;
      }
    }
  
    G4double PRA = PS * Spow;

    if (verboseLevel > 3) {
      cerr << " PS " << PS << " Spow = sqrt(S)*(PR+(1-PQ)*S^4) " << Spow
	     << " PRA = PS*Spow " << PRA <<"\n";
    }

    return std::fabs(PRA);
  }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXInuclParamMomDst_h
