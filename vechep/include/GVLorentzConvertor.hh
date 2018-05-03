//
// $Id: GVLorentzConvertor.hh 2018-05-03 18:34:42Z jlima$
//
// 20180503  Guilherme Lima -- Created, inspired by M.Kelsey's G4LorentzConvertor

#ifndef GXLORENTZ_CONVERTOR_HH
#define GXLORENTZ_CONVERTOR_HH

#include "globals.hh"
#include "LorentzVector.h"
#include "GXTrack.h"
#include "base/vechep/GXThreeVector.h"

//class G4InuclParticle;

using gxbert::GXTrack;
using gxbert::LorentzVector;

template <typename Real_v>
class GVLorentzConvertor {
public:
  GVLorentzConvertor();

  // GVLorentzConvertor(const LorentzVector<Real_v>& bmom, G4double bmass, 
  // 		     const LorentzVector<Real_v>& tmom, G4double tmass);

  GVLorentzConvertor(const GXTrack* bullet, 
		     const GXTrack* target);

  void setVerbose(G4int vb=0) { verboseLevel = vb; }

  void setBullet(const GXTrack* bullet);
  void setTarget(const GXTrack* target);

  void setBullet(const GXTrack& bullet) { setBullet(&bullet); }
  void setTarget(const GXTrack& target) { setTarget(&target); }

  // Use correct four-vectors as input
  void setBullet(const LorentzVector<Real_v>& bmom) {
    bullet_mom = bmom;
    if (verboseLevel > 3) printBullet();
  }

  void setTarget(const LorentzVector<Real_v>& bmom) {
    target_mom = bmom;
    if (verboseLevel > 3) printTarget();
  }

  // These functions "repair" input 4-vectors using specified mass
  void setBullet(const LorentzVector<Real_v>& bmom, G4double bmass) {
    bullet_mom.setVectM(bmom.vect(), bmass);
    if (verboseLevel > 3) printBullet();
  }
  
  void setTarget(const LorentzVector<Real_v>& tmom, G4double tmass) {
    target_mom.setVectM(tmom.vect(), tmass);
    if (verboseLevel > 3) printTarget();
  }

  // Select reference frame for boosts, rotations, etc.
  void toTheCenterOfMass();
  void toTheTargetRestFrame(); 
  void fillKinematics();	// Common calculations after either of above

  LorentzVector backToTheLab(const LorentzVector<Real_v>& mom) const;

  // Four-vectors of bullet and target in last chosen reference frame
  const LorentzVector<Real_v>& getBullet() const { return bullet_mom; }
  const LorentzVector<Real_V>& getTarget() const { return target_mom; }
 
  Real_v getKinEnergyInTheTRS() const;
  Real_v getTotalSCMEnergy() const { return ecm_tot; }
  Real_v getSCMMomentum() const { return scm_momentum.rho(); }
  Real_v getTRSMomentum() const;

  LorentzVector rotate(const LorentzVector& mom) const; 

  LorentzVector rotate(const LorentzVector& mom1,
		       const LorentzVector& mom) const; 

  G4bool reflectionNeeded() const; 

  G4bool trivial() const { return degenerated; }

  // Reporting functions for diagnostics
  void printBullet() const;
  void printTarget() const;

private: 
  static const G4double small;

  G4int verboseLevel;
  LorentzVector bullet_mom;
  LorentzVector target_mom;

  LorentzVector scm_momentum;		// CM momentum relative to target/bullet
  G4ThreeVector   scm_direction;	// Unit vector to reduce repeated calcs

  // Buffer variables for doing ::rotate() calculations
  G4ThreeVector velocity;
  G4double v2;
  G4double ecm_tot;
  G4double valong;
  G4bool degenerated;
};        

#endif // G4LORENTZ_CONVERTOR_HH 
