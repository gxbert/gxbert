//
// $Id: GXLorentzConvertor.hh 2018-05-03 18:34:42Z jlima$
//
// 20180503  Guilherme Lima -- Created, inspired by M.Kelsey's G4LorentzConvertor

#ifndef GXLORENTZ_CONVERTOR_HH
#define GXLORENTZ_CONVERTOR_HH

#include "globals.hh"
#include "LorentzVector.hh"
#include "GXTrack.hh"
#include "GXThreeVector.hh"
#include "GXInuclParticle.hh"
#include "GXHadronicException.hh"
#include "VecCore/VecCore"

using std::cerr;

namespace gxbert {

template <typename Real_v>
class GXLorentzConvertor {

private: 
  int verboseLevel;
  static constexpr double small = 1.0e-10;
  LorentzVector<Real_v> bullet_mom;
  LorentzVector<Real_v> target_mom;

  LorentzVector<Real_v> scm_momentum;	// CM momentum relative to target/bullet
  GXThreeVector<Real_v> scm_direction;	// Unit vector to reduce repeated calcs

  // Buffer variables for doing ::rotate() calculations
  GXThreeVector<Real_v> velocity;
  Real_v v2;
  Real_v ecm_tot;
  Real_v valong;
  vecCore::Mask_v<Real_v> degenerated;

public:
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXLorentzConvertor() {}

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXLorentzConvertor(const LorentzVector<Real_v>& bmom, Real_v bmass,
   		     const LorentzVector<Real_v>& tmom, Real_v tmass);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXLorentzConvertor(const GXTrack &bullet,
  		     const GXTrack &target)
  {
    setBullet(bullet);
    setTarget(target);
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setVerbose(int vb=0) { verboseLevel = vb; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setBullet(Real_v const* px, Real_v const* py, Real_v const* pz, Real_v const* E)
  {
    bullet_mom.x() = px;
    bullet_mom.y() = py;
    bullet_mom.z() = pz;
    bullet_mom.t() = E;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTarget(Real_v const* px, Real_v const* py, Real_v const* pz, Real_v const* E)
  {
    target_mom.x() = px;
    target_mom.y() = py;
    target_mom.z() = pz;
    target_mom.t() = E;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setBullet(const GXTrack &bullet)
  {
    bullet_mom.x() = (Real_v*)bullet.px;
    bullet_mom.y() = (Real_v*)bullet.py;
    bullet_mom.z() = (Real_v*)bullet.pz;
    bullet_mom.t() = (Real_v*)bullet.E;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTarget(const GXTrack &target)
  {
    target_mom.x() = (Real_v*)target.px;
    target_mom.y() = (Real_v*)target.py;
    target_mom.z() = (Real_v*)target.pz;
    target_mom.t() = (Real_v*)target.E;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setBullet(const GXTrack *bullet) { setBullet(*bullet); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTarget(const GXTrack *target) { setTarget(*target); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setBullet(GXInuclParticle<Real_v> const& bullet)
  {
    bullet_mom = bullet.getFourMomentum();
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTarget(GXInuclParticle<Real_v> const& target)
  {
    target_mom = target.getFourMomentum();
  }

  // Use correct four-vectors as input
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setBullet(const LorentzVector<Real_v>& bmom) {
    bullet_mom = bmom;
    if (verboseLevel > 3) printBullet();
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTarget(const LorentzVector<Real_v>& bmom) {
    target_mom = bmom;
    if (verboseLevel > 3) printTarget();
  }

  // These functions "repair" input 4-vectors using specified mass
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setBullet(const LorentzVector<Real_v>& bmom, Real_v bmass) {
    bullet_mom.SetVectMag(bmom.Vect(), bmass);
    if (verboseLevel > 3) printBullet();
  }
  
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setTarget(const LorentzVector<Real_v>& tmom, Real_v tmass) {
    target_mom.SetVectMag(tmom.Vect(), tmass);
    if (verboseLevel > 3) printTarget();
  }

  // Common calculations needed after a boost
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void fillKinematics() {
    ecm_tot = (target_mom+bullet_mom).Mag();

    scm_direction = scm_momentum.Vect().Unit();
    valong = velocity.Dot(scm_direction);

    v2 = velocity.Mag2();

    Real_v pvsq = v2 - valong * valong;  // velocity perp to scm_momentum
    if (verboseLevel > 3) cerr << " pvsq=" << pvsq << "\n";

    degenerated = (pvsq < (Real_v)small);
    if (verboseLevel > 2 && !vecCore::MaskEmpty(degenerated)) 
      cerr << " degenerated case (already along Z): pvsq=" << pvsq <<" small="<< small <<"\n"; 

    if (verboseLevel > 3) {
      cerr << " v2=" << v2 << " valong=" << valong
	     << " valong*valong=" << valong*valong << "\n";
    }
  }

  // Select reference frame for boosts, rotations, etc.
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void toTheCenterOfMass()
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::toTheCenterOfMass" << "\n";

    velocity = (target_mom+bullet_mom).BoostVector();
    if (verboseLevel > 3) cerr << " boost " << velocity << "\n";

    // "SCM" is reverse target momentum in the CM frame
    scm_momentum = target_mom;
    scm_momentum.Boost(-velocity);
    scm_momentum.SetVect(-scm_momentum.Vect());

    if (verboseLevel > 3) cerr << " pscm " << scm_momentum.Vect() << "\n";

    fillKinematics();
  }
  
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void toTheTargetRestFrame()
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::toTheTargetRestFrame" << "\n";

    velocity = target_mom.BoostVector();
    if (verboseLevel > 3) cerr << " boost " << velocity << "\n";

    // "SCM" is bullet momentum in the target's frame
    scm_momentum = bullet_mom;
    scm_momentum.Boost(-velocity);

    if (verboseLevel > 3) cerr << " pseudo-pscm " << scm_momentum.vect() << "\n";

    fillKinematics();
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<Real_v> backToTheLab(const LorentzVector<Real_v>& mom) const
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::backToTheLab" << "\n";

    if (verboseLevel > 3)
      cerr << " at rest: px " << mom.x() << " py " << mom.y() << " pz "
	     << mom.z() << " e " << mom.t() << "\n"
	     << " v2 " << v2 << "\n";

    LorentzVector<Real_v> mom1 = mom;
    mom1.Boost(velocity);

    // use one of the next two options
    auto undo = v2 <= (Real_v)small;
    vecCore::MaskedAssign(mom1.x(), undo, mom.x());
    vecCore::MaskedAssign(mom1.y(), undo, mom.y());
    vecCore::MaskedAssign(mom1.z(), undo, mom.z());
    vecCore::MaskedAssign(mom1.t(), undo, mom.t());

    if (verboseLevel > 3)
      cerr << " at lab: px " << mom1.x() << " py " << mom1.y() << " pz "
	     << mom1.z() << "\n";

    return mom1;
  }

  // Four-vectors of bullet and target in last chosen reference frame
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  const LorentzVector<Real_v>& getBullet() const { return bullet_mom; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  const LorentzVector<Real_v>& getTarget() const { return target_mom; }
 
  // Bullet kinematics in target rest frame (LAB frame, usually)
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Real_v getKinEnergyInTheTRS() const
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::getKinEnergyInTheTRS" << "\n";

    LorentzVector<Real_v> bmom = bullet_mom;
    bmom.Boost(-target_mom.BoostVector());
    return bmom.t()-bmom.Mag();
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Real_v getTotalSCMEnergy() const { return ecm_tot; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Real_v getSCMMomentum() const { return scm_momentum.rho(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Real_v getTRSMomentum() const
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::getTRSMomentum" << "\n";

    LorentzVector<Real_v> bmom = bullet_mom;
    bmom.Boost(-target_mom.boostVector());
    return bmom.rho();
  }

  LorentzVector<Real_v> rotate(const LorentzVector<Real_v>& mom) const
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::rotate(GXLorentzVector)" << "\n";

    if (verboseLevel > 3) {
      cerr << " valong " << valong << " degenerated=" << degenerated << "\n"
	   << " before rotation: px=" << mom.x() << " py=" << mom.y()
	   << " pz=" << mom.z() << "\n";
    }

    LorentzVector<Real_v> mom_rot = mom;
    if (!vecCore::MaskEmpty(degenerated)) {
      if (verboseLevel > 2)
	cerr << " rotating to align with reference z axis " << "\n";

      GXThreeVector<Real_v> vscm = velocity - valong * scm_direction;
      GXThreeVector<Real_v> vxcm = scm_direction.Cross(velocity);

      if (vscm.Mag() > (Real_v)small && vxcm.Mag() > (Real_v)small) {	// Double check
	if (verboseLevel > 3) {
	  cerr << " reference z axis " << scm_direction
	       << " vscm " << vscm << " vxcm " << vxcm << "\n";
	}
      
	mom_rot.SetVect(mom.x()*vscm.Unit() + mom.y()*vxcm.Unit() + mom.z()*scm_direction);
      }
      else {
	if (verboseLevel) 
	  cerr << ">>> GXLorentzVector::rotate zero with !degenerated" << "\n";
      }
    }

    if (verboseLevel > 3) {
      cerr << " after rotation: px " << mom_rot.x() << " py " << mom_rot.y()
	   << " pz " << mom_rot.z() << "\n";
    }

    return mom_rot;
  }

  LorentzVector<Real_v> rotate(const LorentzVector<Real_v>& mom1, const LorentzVector<Real_v>& mom) const
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::rotate(GXLorentzVector,GXLorentzVector)"
	   << "\n";

    if (verboseLevel > 3) {
      cerr << " before rotation: px " << mom.x() << " py " << mom.y()
	   << " pz " << mom.z() << "\n";
    }

    GXThreeVector<Real_v> mom1_dir = mom1.Vect().Unit();
    Real_v pv = velocity.Dot(mom1_dir);

    Real_v vperp = v2 - pv*pv;		// velocity perpendicular to mom1
    if (verboseLevel > 3) {
      cerr << " vperp " << vperp << " small? " << (vperp <= small) << "\n";
    }

    LorentzVector<Real_v> mom_rot = mom;

    if (vperp > (Real_v)small) {
      if (verboseLevel > 2)
	cerr << " rotating to align with first z axis " << "\n";

      GXThreeVector<Real_v> vmom1 = velocity - pv*mom1_dir;
      GXThreeVector<Real_v> vxm1  = mom1_dir.Cross(velocity);

      const Real_v small2 = small*small;
      if (vmom1.Mag2() > small2 && vxm1.Mag2() > small2) {	// Double check
	if (verboseLevel > 3) {
	  cerr << " first z axis " << mom1_dir << "\n"
	       << " vmom1 " << vmom1 << " vxm1 " << vxm1 << "\n";
	}
      
	mom_rot.SetVect(mom.x()*vmom1.Unit() + mom.y()*vxm1.Unit() +
			mom.z()*mom1_dir );
      }
      else {
	if (verboseLevel)
	  cerr << ">>> GXLorentzVector::rotate zero with !degenerated" << "\n";
      }
    }

    if (verboseLevel > 3) {
      cerr << " after rotation: px " << mom_rot.x() << " py " << mom_rot.y()
	   << " pz " << mom_rot.z() << "\n";
    }

    return mom_rot;
  }

  vecCore::Mask_v<Real_v> reflectionNeeded() const
  {
    if (verboseLevel > 2)
      cerr << " >>> GXLorentzConvertor::reflectionNeeded (query)" << "\n";

    if (verboseLevel > 3) {
      cerr << " v2 = " << v2 << " SCM z = " << scm_momentum.z()
	   << " degenerated? " << degenerated << "\n";
    }

    if (v2 < (Real_v)small && !vecCore::MaskFull(degenerated))
      throw GXHadronicException(__FILE__, __LINE__, "GXLorentzConvertor::reflectionNeeded - return value undefined");

    if (verboseLevel > 2) {
      cerr << " reflection across XY is"
	   << ((v2>=(Real_v)small && (!vecCore::MaskFull(degenerated) || scm_momentum.z()<0.0))?"":" NOT")
	   << " needed" << "\n";
    }

    return (v2>=(Real_v)small && (!vecCore::MaskFull(degenerated) || scm_momentum.z()<0.0));
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  vecCore::Mask_v<Real_v> trivial() const { return degenerated; }

  // Reporting functions for diagnostics
  void printBullet() const
  {
    cerr << " GXLC bullet: px " << bullet_mom.x() << " py " << bullet_mom.y()
	 << " pz " << bullet_mom.z() << " e " << bullet_mom.t()
	 << " mass " << bullet_mom.Mag() << "\n";
  }

  void printTarget() const
  {
    cerr << " GXLC target: px " << target_mom.x() << " py " << target_mom.y()
	 << " pz " << target_mom.z() << " e " << target_mom.t()
	 << " mass " << target_mom.Mag() << "\n";
  }

};        

} // end namespace gxbert

#endif // GXLORENTZ_CONVERTOR_HH 
