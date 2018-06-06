//
// @File: GXInuclElementaryParticle.h
//
// 20180528  Guilherme Lima -- Created, based on M.Kelsey's G4InuclElementaryParticle

#ifndef GXBERT_GXInuclElementaryParticle_H
#define GXBERT_GXInuclElementaryParticle_H

//#include "VecCore/VecCore"
//#include "IntFor.hh"
#include "ApproxEqual.hh"
#include "GXThreeVector.hh"
#include "LorentzVector.hh"
#include "G4NucleiModel.hh"
#include "GXParticleDefinition.hh"
#include "G4InuclParticleNames.hh"

//#include "VecRng/MRG32k3a.h"

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
#include <iostream>
#include <strstream>

using namespace G4InuclParticleNames;

namespace gxbert {

inline namespace GXBERT_IMPL_NAMESPACE {

template <typename T>
class GXInuclElementaryParticle {

  using Bool_v = vecCore::Mask_v<T>;

public:
  // used to indicate model that created instance of G4InuclParticle
  // 0 default
  // 1 bullet
  // 2 target
  // 3 G4ElementaryParticleCollider
  // 4 G4IntraNucleiCascader
  // 5 G4NonEquilibriumEvaporator
  // 6 G4EquilibriumEvaporator
  // 7 G4Fissioner
  // 8 G4BigBanger
  // 9 G4PreCompound
  // 10 G4CascadeCoalescence
  enum Model { DefaultModel, bullet, target, EPCollider, INCascader,
	       NonEquilib, Equilib, Fissioner, BigBanger, PreCompound,
	       Coalescence };

private:
  GXThreeVector<T> fDir;
  T fkinEnergy;
  Model fModelID;
  Index_v<T> iType;
  //IntFor<T> iType;
  //vecRng::MRG32k3a<ScalarBackend> *fRNG;

public:
  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  static GXParticleDefinition const* getDefinition(int ityp)
  {
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
      std::cerr << "GXInuclElementaryParticle::getDefinition(): unknown particle type "
		<< ityp << G4endl;
    }

    return 0;
  }

  /// Constructors
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle()
    : fDir(0.,0.,1.)
    , fkinEnergy(0.)
    , fModelID(DefaultModel)
    , iType(-1)
  { }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(Index_v<T> ityp, Model model = DefaultModel)
    : fDir(0.,0.,1.)
    , fkinEnergy(0.)
    , fModelID(model)
    , iType(ityp)
  { }

  /*VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const GXDynamicParticle& dynPart,
			    Model model = DefaultModel)
    : fDir(dynPart.fDir)
    , fkinEnergy(dynPart.fkinEnergy)
    , fModelID(model)
    , iType(dynPart.getType())
    { }*/

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const LorentzVector<T>& mom,
			    Index_v<T> const& ityp, Model model = DefaultModel)
    : fDir(mom.Vect().Unit())
    , fkinEnergy(mom.t() - mom.Mag())
    , fModelID(model)
    , iType(ityp)
  { }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(T ekin, Index_v<T> ityp,
			    Model model = DefaultModel) 
    : fDir(0., 0., 1.)
    , fkinEnergy(ekin)
    , fModelID(model)
    , iType(ityp)
  { }

  // WARNING:  This may create a particle without a valid type code!
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const LorentzVector<T>& mom,
			    const GXParticleDefinition* pd,
			    Model model=DefaultModel)
    : fDir(mom.Vect().Unit())
    , fkinEnergy(mom.t() - mom.Mag())
    , fModelID(model)
    , iType(pd->GetParticleType())
  { }

  // Copy constructor
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const GXInuclElementaryParticle<T>& right)
    : fDir(right.fDir)
    , fkinEnergy(right.fkinEnergy)
    , fModelID(right.fModelID)
    , iType(right.iType)
  { }

  // assignment operator
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  G4InuclElementaryParticle& operator=(const GXInuclElementaryParticle<T>& right)
  {
    fDir.Set(right.x(), right.y(), right.z());
    fModelID = right.getModel();
    fkinEnergy = right.fkinEnergy;
    iType = right.iType;
  }

  // comparison
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<T> operator==(const GXInuclElementaryParticle<T>& right)
  {
    return (fDir == right.fDir && ApproxEqual(fkinEnergy, right.fkinEnergy) && fModelID == right.fModelID && iType == right.iType);
  }
  
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<T> operator!=(const GXInuclElementaryParticle<T>& right)
  {
    return ! (this->operator==(right));
  }
  
  ~GXInuclElementaryParticle() {}

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void fill(const LorentzVector<T>& mom, Index_v<T> const& itype, Model model = DefaultModel)
  {
    fDir = mom.Vect().Unit();
    fkinEnergy = mom.t() - mom.Mag();
    fModelID = model;
    iType = itype;
  }

  // Ensure that type code refers to a known particle
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  static Bool_v valid(Index_v<T> const& ityp) { return ityp != 0; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v valid() const { return valid(type()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v isPhoton() const { return isOfType(photon); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v isMuon() const { return isOfType(mup) | isOfType(mum); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v isElectron() const { return Bool_v(G4InuclParticleNames::isElectron(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v isNeutrino() const { return Bool_v(G4InuclParticleNames::isNeutrino(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v pion() const { return isOfType(pionPlus) | isOfType(pionMinus) | isOfType(pionZero); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v nucleon() const { return isOfType(proton) | isOfType(neutron); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v antinucleon() const { return Bool_v(G4InuclParticleNames::antinucleon(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Index_v<T> baryon() const
  {
    // Can use as a bool (!=0 ==> true)
    Index_v<T> result(0);
    for (size_t i = 0; i < VectorSize<T>(); ++i)
      Set(result, i, getDefinition(Get(iType, i))->GetBaryonNumber());
    return result;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v antibaryon() const { return baryon() < 0; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v hyperon() const { return (baryon() && getStrangeness()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v quasi_deutron() const {
    return (iType > Index_v<T>(100));
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Index_v<T> getStrangeness() const { return getStrangeness(type()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void print(std::ostream& os) const
  {
    os << "\n" << " Particle: " << getDefinition(type())->GetParticleName()
       << " type=" << type() << " mass=" << getParticleMass()
       << " ekin=" << getKineticEnergy();
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  static Index_v<T> getStrangeness(Index_v<T> type) {
    Index_v<T> result(0);
    for (size_t i = 0; i < VectorSize<T>(); ++i)
      Set(result, i, getDefinition(Get(type, i))->GetStrangeness());
    return result;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  static T getParticleMass(Index_v<T> types)
  {
    T result;
    for (size_t i = 0; i < VectorSize<T>(); ++i) {
      auto pd = getDefinition(Get(types, i));
      Set(result, i, (pd ? pd->GetPDGMass() : 0.0));
    }
    return result;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setKineticEnergy(T const& ekin) { fkinEnergy = ekin; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T getKineticEnergy() const { return fkinEnergy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T getTotalEnergy() const { return fkinEnergy + getParticleMass(iType); }

  // VECCORE_ATT_HOST_DEVICE
  // VECCORE_FORCE_INLINE
  // void setType(IntFor<T> const& type) { iType = type; }

  // VECCORE_ATT_HOST_DEVICE
  // VECCORE_FORCE_INLINE
  // IntFor<T> type() { return iType; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setType(Index<T> const& type) { iType = type; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Index<T> type() const { return iType; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v isOfType(int ityp) const
  {
    Index_v<T> ref( (Index<T>)ityp );
    return iType == ref;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setMomentumDirection(GXThreeVector<T> const& momentum)
  {
    fDir = momentum;
    assert(momentum.Mag2() == 1.0);
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setModel(Index<T> const& model)
  {
    fModelID = model;
  }

  static int type(const GXParticleDefinition *pd) {
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
      //return ((inuclRndm() > 0.5) ? kaonZero : kaonZeroBar);
      //return (fRNG->Uniform<ScalarBackend>() > 0.5 ? kaonZero : kaonZeroBar);
      return kaonZero; // ??? TODO: fix this
    }

    return 0;	// Unknown objects return zero (e.g., nuclei)
  }

  // Overwrite data structure (avoids creating/copying temporaries)

  void fill(LorentzVector<T> const& mom, Index_v<T> const& ityp, Model const& model)
  {
    setType(ityp);
    setMomentum(mom);
    setModel(model);
  }

  void fill(T const& ekin, Index_v<T>& ityp, Model const& model)
  {
    setType(ityp);
    setKineticEnergy(ekin);
    setModel(model);
  }

  void fill(LorentzVector<T> const& mom, const GXParticleDefinition* pd, Model const& model)
  {
    //setDefinition(pd);
    setMomentum(mom);
    setModel(model);
  }

};

template <typename T>
std::ostream &operator<<(std::ostream &os, GXInuclElementaryParticle<T> const& trk)
{
  auto vsize = vecCore::VectorSize<T>();
  std::strstream names;
  std::strstream masses;
  if (vsize > 1) {
    auto ityp = trk.type();
    for (int i = 0; i < vsize; ++i) {
      if (i > 0) {
	names <<"; ";
	masses <<"; ";
      }
      GXParticleDefinition const* pd = trk.getDefinition(Get(ityp, i));
      names << pd->GetParticleName();
      masses << pd->GetPDGMass();
    }
    os << " Particles=["<< names.str() <<"]"
       << " masses=[" << masses.str() <<"]"
       << " types=" << trk.type()
       << " ekin=" << trk.getKineticEnergy();
  }
  return os;
}

template <>
std::ostream &operator<<(std::ostream &os, GXInuclElementaryParticle<double> const& trk)
{
  GXParticleDefinition const* pd = trk.getDefinition(trk.type());
  os << " Particle: " << pd->GetParticleName()
     << " mass=" << pd->GetPDGMass()
     << " type=" << trk.type()
     << " ekin=" << trk.getKineticEnergy();
  return os;
}

} // end GXBERT_IMPL_NAMESPACE
} // end namespace gxbert
#endif // GXBERT_GXInuclElementaryParticle_H
