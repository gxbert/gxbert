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
  //GXParticleDefinition *fParticleDef;
  Index_v<T> iType;
  //IntFor<T> iType;

public:
  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  static GXParticleDefinition const* makeDefinition(int type);

  /// Constructors
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle()
    : fDir(0.,0.,1.)
    , fkinEnergy(0.)
    , fModelID(DefaultModel)
      //, fParticleDef(NULL)
    , iType(-1)
  { }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(int ityp, Model model = DefaultModel) 
    : fDir(0.,0.,1.)
    , fkinEnergy(0.)
    , fModelID(model)
      //, fParticleDef(GXInuclElementaryParticle::makeDefinition(ityp))
    , iType(ityp)
  { }

  /*VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const GXDynamicParticle& dynPart,
			    Model model = DefaultModel)
    : fDir(dynPart.fDir)
    , fkinEnergy(dynPart.fkinEnergy)
    , fModelID(model)
    //, fParticleDef(dynPart.fParticleDef)
    , iType(dynPart.getType())
    { }*/

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const LorentzVector<T>& mom,
			    int ityp, Model model = DefaultModel)
    : fDir(mom.Vect().Unit())
    , fkinEnergy(mom.t() - mom.Mag())
    , fModelID(model)
      //, fParticleDef(makeDefinition(ityp))
    , iType(ityp)
  { }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(double ekin, int ityp,
			    Model model = DefaultModel) 
    : fDir(0., 0., 1.)
    , fkinEnergy(ekin)
    , fModelID(model)
      //, fParticleDef(makeDefinition(ityp))
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
    // , fParticleDef(pd)
    , iType(pd->GetParticleType())
  { }

  // Copy constructor
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclElementaryParticle(const GXInuclElementaryParticle& right)
    : fDir(right.fDir)
    , fkinEnergy(right.fkinEnergy)
    , fModelID(right.fModelID)
    // , fParticleDef(other.getParticleDefinition())
    , iType(right.iType)
  { }

  // assignment operator
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  G4InuclElementaryParticle& operator=(const GXInuclElementaryParticle& right)
  {
    fDir.Set(right.x(), right.y(), right.z());
    fModelID = right.getModel();
    fkinEnergy = right.fkinEnergy;
    //fParticleDef = right.getParticleDefinition();
    iType = right.iType;
  }

  // comparison
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<T> operator==(const GXInuclElementaryParticle& right)
  {
    return (fDir == right.fDir && ApproxEqual(fkinEnergy, right.fkinEnergy) && fModelID == right.fModelID && iType == right.iType);
  }
  
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<T> operator!=(const GXInuclElementaryParticle& right)
  {
    return ! (this->operator==(right));
  }
  
  ~GXInuclElementaryParticle() {}

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void fill(const LorentzVector<T>& mom, const Int_v itype, Model model = DefaultModel)
  {
    fDir = mom.Vect().Unit();
    fkinEnergy = mom.t() - mom.Mag();
    fModelID = model;
    iType = itype;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v isPhoton() const { return Bool_v(G4InuclParticleNames::isPhoton(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v isMuon() const { return Bool_v(G4InuclParticleNames::isMuon(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v isElectron() const { return Bool_v(G4InuclParticleNames::isElectron(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v isNeutrino() const { return Bool_v(G4InuclParticleNames::isNeutrino(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v pion() const { return Bool_v(G4InuclParticleNames::pion(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v nucleon() const { return Bool_v(G4InuclParticleNames::nucleon(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
    Bool_v antinucleon() const { return Bool_v(G4InuclParticleNames::antinucleon(type())); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setKineticEnergy(T const& ekin) { fkinEnergy = ekin; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T getKineticEnergy() const { return fkinEnergy; }

  // VECCORE_ATT_HOST_DEVICE
  // VECCORE_FORCE_INLINE
  // void setType(IntFor<T> const& type) { /*fParticleDef = makeDefinition(type);*/ iType = type; }

  // VECCORE_ATT_HOST_DEVICE
  // VECCORE_FORCE_INLINE
  // IntFor<T> type() { return iType; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setType(Index<T> const& type) { /*fParticleDef = makeDefinition(type);*/ iType = type; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Index<T> type() { return iType; }

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
};

} // end GXBERT_IMPL_NAMESPACE
} // end namespace gxbert
#endif // GXBERT_GXInuclElementaryParticle_H
