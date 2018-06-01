//
// @File: GXInuclElementaryParticle.cc
//
// 20180530 Guilherme Lima -- Created, based on M.Kelsey's G4InuclElementaryParticle class

#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"
using namespace G4InuclSpecialFunctions;

#include "GXSystemOfUnits.hh"
#include "GXParticleDefinition.hh"
//#include "G4ParticleTypes.hh"

#include "GXProton.hh"
#include "GXNeutron.hh"
#include "GXPionPlus.hh"
#include "GXPionMinus.hh"
#include "GXPionZero.hh"
#include "GXGamma.hh"
#include "GXKaonPlus.hh"
#include "GXKaonMinus.hh"
#include "GXKaonZero.hh"
#include "GXKaonZeroLong.hh"
#include "GXKaonZeroShort.hh"
#include "GXAntiKaonZero.hh"
#include "GXLambda.hh"
#include "GXSigmaPlus.hh"
#include "GXSigmaZero.hh"
#include "GXSigmaMinus.hh"
#include "GXXiZero.hh"
#include "GXXiMinus.hh"
#include "GXOmegaMinus.hh"
#include "GXDeuteron.hh"
#include "GXTriton.hh"
#include "GXHe3.hh"
#include "GXAlpha.hh"

#include "G4Diproton.hh"
#include "G4UnboundPN.hh"
#include "G4Dineutron.hh"
//#include "Randomize.hh"

#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

const GXParticleDefinition*
G4InuclElementaryParticle::makeDefinition(G4int ityp) {
  switch(ityp) {
  case proton:      return GXProton::Definition(); break;
  case neutron:     return GXNeutron::Definition(); break;
  case pionPlus:    return GXPionPlus::Definition(); break;
  case pionMinus:   return GXPionMinus::Definition(); break;
  case pionZero:    return GXPionZero::Definition(); break;
  case photon:      return GXGamma::Definition(); break;
  case kaonPlus:    return GXKaonPlus::Definition(); break;
  case kaonMinus:   return GXKaonMinus::Definition(); break;
  case kaonZero:    return GXKaonZero::Definition(); break;
  case kaonZeroBar: return GXAntiKaonZero::Definition(); break;
  case lambda:      return GXLambda::Definition(); break;
  case sigmaPlus:   return GXSigmaPlus::Definition(); break;
  case sigmaZero:   return GXSigmaZero::Definition(); break;
  case sigmaMinus:  return GXSigmaMinus::Definition(); break;
  case xiZero:      return GXXiZero::Definition(); break;
  case xiMinus:     return GXXiMinus::Definition(); break;
  case omegaMinus:  return GXOmegaMinus::Definition(); break;
    // NOTE:  The four light nuclei "particles" are actually G4Ions
  case deuteron:    return GXDeuteron::Definition(); break;
  case triton:      return GXTriton::Definition(); break;
  case He3:	    return GXHe3::Definition(); break;
  case alpha:	    return GXAlpha::Definition(); break;
/*
  case antiProton:  return G4AntiProton::Definition(); break;
  case antiNeutron: return G4AntiNeutron::Definition(); break;
    // NOTE:  The the four light antinuclei "particles" are actually G4Ions
  case antiDeuteron: return G4AntiDeuteron::Definition(); break;
  case antiTriton:  return G4AntiTriton::Definition(); break;
  case antiHe3:     return G4AntiHe3::Definition(); break;
  case antiAlpha:   return G4AntiAlpha::Definition(); break;
    // NOTE:  The three unbound dibaryons are local Bertini classes
*/
  case diproton:    return G4Diproton::Definition(); break;
  case unboundPN:   return G4UnboundPN::Definition(); break;
  case dineutron:   return G4Dineutron::Definition(); break;
/*
    // Leptons are included for muon capture and future tau/neutrino physics
  case electron:    return G4Electron::Definition(); break;
  case positron:    return G4Positron::Definition(); break;
  case electronNu:  return G4NeutrinoE::Definition(); break;
  case antiElectronNu: return G4AntiNeutrinoE::Definition(); break;
  case muonMinus:   return G4MuonMinus::Definition(); break;
  case muonPlus:    return G4MuonPlus::Definition(); break;
  case muonNu:      return G4NeutrinoMu::Definition(); break;
  case antiMuonNu:  return G4AntiNeutrinoMu::Definition(); break;
  case tauMinus:    return G4TauMinus::Definition(); break;
  case tauPlus:     return G4TauPlus::Definition(); break;
  case tauNu:       return G4NeutrinoTau::Definition(); break;
  case antiTauNu:   return G4AntiNeutrinoTau::Definition(); break;
*/
  default:
    G4cerr << "G4InuclElementaryParticle::makeDefinition: unknown particle type "
           << ityp << G4endl;
  }
  
  return 0;
}

// This is the inverse mapping to makeDefinition above

G4int G4InuclElementaryParticle::type(const GXParticleDefinition *pd) {
  if (pd == 0) return 0;
  if (pd == GXProton::Definition())       return proton;
  if (pd == GXNeutron::Definition())      return neutron;
  if (pd == GXPionPlus::Definition())     return pionPlus;
  if (pd == GXPionMinus::Definition())    return pionMinus;
  if (pd == GXPionZero::Definition())     return pionZero;
  if (pd == GXGamma::Definition())        return photon;
  if (pd == GXKaonPlus::Definition())     return kaonPlus;
  if (pd == GXKaonMinus::Definition())    return kaonMinus;
  if (pd == GXKaonZero::Definition())     return kaonZero;
  if (pd == GXAntiKaonZero::Definition()) return kaonZeroBar;
  if (pd == GXLambda::Definition())       return lambda;
  if (pd == GXSigmaPlus::Definition())    return sigmaPlus;
  if (pd == GXSigmaZero::Definition())    return sigmaZero;
  if (pd == GXSigmaMinus::Definition())   return sigmaMinus;
  if (pd == GXXiZero::Definition())       return xiZero;
  if (pd == GXXiMinus::Definition())      return xiMinus;
  if (pd == GXOmegaMinus::Definition())   return omegaMinus;
  // NOTE:  The four light nuclei "particles" are actually G4Ions
  if (pd == GXDeuteron::Definition())     return deuteron;
  if (pd == GXTriton::Definition())       return triton;
  if (pd == GXHe3::Definition())          return He3;
  if (pd == GXAlpha::Definition())        return alpha;
/*
  if (pd == G4AntiProton::Definition())   return antiProton;
  if (pd == G4AntiNeutron::Definition())  return antiNeutron;
  // NOTE:  The the four light antinuclei "particles" are actually G4Ions
  if (pd == G4AntiDeuteron::Definition()) return antiDeuteron;
  if (pd == G4AntiTriton::Definition())   return antiTriton;
  if (pd == G4AntiHe3::Definition())      return antiHe3;
  if (pd == G4AntiAlpha::Definition())    return antiAlpha;
*/
  // NOTE:  The three unbound dibaryons are local Bertini classes
  if (pd == G4Diproton::Definition())     return diproton;
  if (pd == G4UnboundPN::Definition())    return unboundPN;
  if (pd == G4Dineutron::Definition())    return dineutron;

/*
  if (pd == G4Electron::Definition())     return electron;
  if (pd == G4Positron::Definition())     return positron;
  if (pd == G4NeutrinoE::Definition())    return electronNu;
  if (pd == G4AntiNeutrinoE::Definition()) return antiElectronNu;
  if (pd == G4MuonMinus::Definition())    return muonMinus;
  if (pd == G4MuonPlus::Definition())     return muonPlus;
  if (pd == G4NeutrinoMu::Definition())   return muonNu;
  if (pd == G4AntiNeutrinoMu::Definition()) return antiMuonNu;
  if (pd == G4TauMinus::Definition())     return tauMinus;
  if (pd == G4TauPlus::Definition())      return tauPlus;
  if (pd == G4NeutrinoTau::Definition())  return tauNu;
  if (pd == G4AntiNeutrinoTau::Definition()) return antiTauNu;

  // Weak neutral kaons must be mixed back to strong (strangeness states)
*/
  if (pd==GXKaonZeroShort::Definition() || pd==GXKaonZeroLong::Definition()) {
    return ((inuclRndm() > 0.5) ? kaonZero : kaonZeroBar);
  }

  return 0;	// Unknown objects return zero (e.g., nuclei)
}

void G4InuclElementaryParticle::setType(G4int ityp) {
  setDefinition(makeDefinition(ityp));
}


// Overwrite data structure (avoids creating/copying temporaries)

void G4InuclElementaryParticle::fill(const G4LorentzVector& mom, G4int ityp,
				     G4InuclParticle::Model model) {
  setType(ityp);
  setMomentum(mom);
  setModel(model);
}

void G4InuclElementaryParticle::fill(G4double ekin, G4int ityp,
				     G4InuclParticle::Model model) {
  setType(ityp);
  setKineticEnergy(ekin);
  setModel(model);
}

void G4InuclElementaryParticle::fill(const G4LorentzVector& mom,
				     const GXParticleDefinition* pd,
				     G4InuclParticle::Model model) {
  setDefinition(pd);
  setMomentum(mom);
  setModel(model);
}


// Assignment operator for use with std::sort()
G4InuclElementaryParticle& 
G4InuclElementaryParticle::operator=(const G4InuclElementaryParticle& right) {
  G4InuclParticle::operator=(right);
  return *this;
}


G4int G4InuclElementaryParticle::getStrangeness(G4int ityp) {
  const GXParticleDefinition* pd = makeDefinition(ityp);
  return pd ? (pd->GetQuarkContent(3) - pd->GetAntiQuarkContent(3)) : 0;
}

G4double G4InuclElementaryParticle::getParticleMass(G4int ityp) {
  const GXParticleDefinition* pd = makeDefinition(ityp);
  return pd ? pd->GetPDGMass()*MeV/GeV : 0.0;	// From G4 to Bertini units
}


// Print particle parameters

void G4InuclElementaryParticle::print(std::ostream& os) const {
  G4InuclParticle::print(os);
  os << G4endl << " Particle: " << getDefinition()->GetParticleName() 
     << " type " << type() << " mass " << getMass()
     << " ekin " << getKineticEnergy(); 
}

