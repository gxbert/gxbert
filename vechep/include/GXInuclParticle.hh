//
// @File: GXInuclParticle.h
//
// 20180528  Guilherme Lima -- Created, based on M.Kelsey's G4InuclParticle

#ifndef GXBERT_GXInuclParticle_H
#define GXBERT_GXInuclParticle_H

#include "ApproxEqual.hh"
#include "GXThreeVector.hh"
#include "LorentzVector.hh"
//#include "G4NucleiModel.hh"
#include "GXParticleDefinition.hh"
#include "G4InuclParticleNames.hh"

//#include "VecRng/MRG32k3a.h"

//#include <iostream>

using namespace G4InuclParticleNames;

namespace gxbert {

inline namespace GXBERT_IMPL_NAMESPACE {

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

template <typename T>
class GXInuclParticle {

using Bool_v = vecCore::Mask_v<T>;

protected:
  GXThreeVector<T> fDir;
  T fkinEnergy;
  T fMass;
  Model fModelID;

public:

  /// Constructors

  // Note: a default constructor is implicit here
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle(Model model = DefaultModel)
    : fDir(0.,0.,1.)
    , fkinEnergy(0.)
    , fModelID(model)
  { }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle(const LorentzVector<T>& mom, Model model = DefaultModel)
    : fDir(mom.Vect().Unit())
    , fkinEnergy(mom.t() - mom.Mag())
    , fModelID(model)
  { }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle(T ekin, Index_v<T> ityp, Model model = DefaultModel) 
    : fDir(0., 0., 1.)
    , fkinEnergy(ekin)
    , fModelID(model)
  { }

  // Copy constructor
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle(const GXInuclParticle<T>& right)
    : fDir(right.fDir)
    , fkinEnergy(right.fkinEnergy)
    , fModelID(right.fModelID)
  { }

  // assignment operation is fully defined in derived classes
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle& operator=(const GXInuclParticle<T>& right) = delete;

  // comparison operations are fully defined in derived classes
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<T> operator==(const GXInuclParticle<T>& right)
  {
    assert(false);
  }
  
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<T> operator!=(const GXInuclParticle<T>& right)
  {
    assert(false);
  }
  
  ~GXInuclParticle() {}


  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setKineticEnergy(T const& ekin) { fkinEnergy = ekin; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T getKineticEnergy() const { return fkinEnergy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTotalEnergy(T const& etot) { fkinEnergy = etot - fMass; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T getTotalEnergy() const { return fkinEnergy + fMass; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setMomentumDirection(GXThreeVector<T> const& momentum)
  {
    fDir = momentum;
    assert( vecCore::MaskFull( ApproxEqual( momentum.Mag2(), T(1.0)) ) );
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setMomentum(LorentzVector<T> const& lorvec)
  {
    fDir = lorvec.Vect().Unit();
    setTotalEnergy( lorvec.E() );
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector<T> const& getMomentumDirection() const
  {
    return fDir;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T getMomModule() const
  {
    return vecCore::math::Sqrt(fkinEnergy * ( fkinEnergy + T(2.0) * fMass));
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector<T> getMomentum() const
  {
    return getMomModule() * fDir;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T> getFourMomentum() const
  {
    LorentzVector<T> temp;
    T momValue = getMomModule();
    temp.Set( momValue * fDir.x(), momValue * fDir.y(), momValue * fDir.z(), fkinEnergy + fMass );
    return temp;
  }

protected:
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector<T>& getMomentumDirection()
  {
    return fDir;
  }

public:
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Model getModel() const
  {
    return fModelID;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setModel(Model const& model)
  {
    fModelID = model;
  }

  /// retrieves PDG masses based on particle types provided as argument
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  virtual void setParticleMass(Index_v<T> const& itype) = 0;

  /// Mass setter
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setMass(T const& mass) { fMass = mass; }

  /// Mass getters
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& getParticleMass() const { return fMass; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& mass() const { return fMass; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  virtual void print(std::ostream& os) const = 0;
};

} // end GXBERT_IMPL_NAMESPACE
} // end namespace gxbert
#endif // GXBERT_GXInuclParticle_H
