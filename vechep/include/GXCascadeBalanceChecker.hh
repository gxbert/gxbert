//
// @File: GXCascadeBalanceChecker.hh
//
// 20181004  Guilherme Lima -- Created, based on M.Kelsey's class G4CascadeCheckBalance

#ifndef GXCascadeBalanceChecker_hh
#define GXCascadeBalanceChecker_hh

//#include "GXVCascadeCollider.hh"
//#include "GXCollisionOutput.hh"
//#include "LorentzVector.hh"
#include <cmath>
#include <vector>

////class G4CascadParticle;
//#include "GXInuclParticle.hh"
//#include "GXInuclElementaryParticle.hh"
//#include "G4InuclNuclei.hh"

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

template <typename T>
class GXCascadeBalanceChecker {

private: // typenames

  using Bool_v = Mask_v<T>;

  constexpr static size_t fvsize = VectorSize<T>();
  using Int_v = typename vecCore::backend::VcSimdArray<fvsize>::Int_v;

public:
  static const T tolerance;	// Don't do floating zero!

  explicit GXCascadeBalanceChecker(const std::string& name = "GXCascadeBalanceChecker")
    : verboseLevel(0)
    , fName(name)
  { }

  GXCascadeBalanceChecker(T relative, T absolute, const std::string& owner = "GXCascadeBalanceChecker");

  virtual ~GXCascadeBalanceChecker() { }

  void setVerboseLevel(int verbose) { verboseLevel = verbose; }

  void setLimits(T relative, T absolute) {
    setRelativeLimit(relative);
    setAbsoluteLimit(absolute);
  }

  void setRelativeLimit(T limit) { relativeLimit = limit; }
  void setAbsoluteLimit(T limit) { absoluteLimit = limit; }

  void check(GXInuclParticle<T> const* bullet,
	     GXInuclParticle<T> const* target,
	     GXCollisionOutput<T> const& output);

  void check(GXInuclParticle<T> const* bullet,
	     GXInuclParticle<T> const* target,
	     std::vector<GXInuclElementaryParticle<T>> const& daughters);

  /*
  // This is for use with G4VCascadeDeexcitation modules
  void collide(const GXFragment& fragment, GXCollisionOutput& output);

  // This is for use with GXEPCollider internal checks
  void collide(GXInuclParticle* bullet, GXInuclParticle* target,
	       const std::vector<GXInuclElementaryParticle>& particles);

  // This is for use with GXNucleiModel internal checks
  void collide(GXInuclParticle* bullet, GXInuclParticle* target,
	       const std::vector<GXCascadParticle>& particles);

  // This is for use with GXIntraNucleiCascader
  void collide(GXInuclParticle* bullet, GXInuclParticle* target,
	       GXCollisionOutput& output,
	       const std::vector<GXCascadParticle>& cparticles);

  // This is for use with G4BigBanger internal checks
  void collide(const GXFragment& target,
	       const std::vector<GXInuclElementaryParticle>& particles);

  // This is for use with G4Fissioner internal checks
  void collide(const GXFragment& target,
	       const std::vector<G4InuclNuclei>& fragments);
  */

  // Checks on conservation laws (kinematics, baryon number, charge, hyperons)
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v energyOkay() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v ekinOkay() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v momentumOkay() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v baryonOkay() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v chargeOkay() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v strangeOkay() const;

  // Global check, used by GXCascadeInterface validation loop
  // NOTE:  Strangeness is not required to be conserved in final state
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v okay() const
  {
    return (energyOkay() && momentumOkay() && baryonOkay() && chargeOkay());
  }

  // Calculations of conserved quantities from initial and final state
  // FIXME:  Relative comparisons don't work for zero!
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T deltaE() const { return (final.E() - initial.E()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T relativeE() const {
    T result(0.);
    Mask_v<T> done(math::Abs(deltaE()) < tolerance);
    MaskedAssign(result, !done, Blend((initial.E() < tolerance), T(1.), deltaE()/initial.E() ));
    return result;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T deltaKE() const { return (ekin(final) - ekin(initial)); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T relativeKE() const {
    return ( (math::Abs(deltaKE()) < tolerance) ? 0. : 
	     (ekin(initial) < tolerance) ? 1. : deltaKE()/ekin(initial) );
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T deltaP() const { return deltaLV().Perp(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T relativeP() const {
    return Blend( (math::Abs(deltaP()) < tolerance), T(0.),
		  Blend((initial.Perp2() < tolerance*tolerance), T(1.), deltaP()/initial.Perp() ) );
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T> deltaLV() const { return final - initial; }

  // Baryon number, charge, S are discrete; no bounds and no "relative" scale
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Int_v deltaB() const { return (finalBaryon - initialBaryon); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Int_v deltaQ() const { return (finalCharge - initialCharge); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Int_v deltaS() const { return (finalStrange- initialStrange); }

protected:
  // Utility function for kinetic energy
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T ekin(const LorentzVector<T>& p) const { return (p.e() - p.m()); }

private:
  int verboseLevel;
  std::string fName;

  T relativeLimit;	// Fractional bound on conservation
  T absoluteLimit;	// Absolute (GeV) bound on conservation

  LorentzVector<T> initial;	// Four-vectors for computing violations
  LorentzVector<T> final;

  Int_v initialBaryon;		// Total baryon number
  Int_v finalBaryon;

  Int_v initialCharge;		// Total charge
  Int_v finalCharge;

  Int_v initialStrange;		// Total strangeness (s-quark content)
  Int_v finalStrange;

  GXCollisionOutput<T> tempOutput;		// Buffer for direct-list interfaces

private:
  // Copying of modules is forbidden
  GXCascadeBalanceChecker(const GXCascadeBalanceChecker&);
  GXCascadeBalanceChecker& operator=(const GXCascadeBalanceChecker&);
};

  template <typename T>
  const T GXCascadeBalanceChecker<T>::tolerance(1e-6);	// How small is zero?

  template <typename T>
  void GXCascadeBalanceChecker<T>::
  check(GXInuclParticle<T> const* bullet, GXInuclParticle<T> const* target,
	GXCollisionOutput<T> const& output)
  {
    if (verboseLevel>1) std::cerr << " >>> GXCascadeBalanceChecker(" << fName << ")::check()...\n";

    initial *= 0.;	// Fast reset; some colliders only have one pointer
    if (bullet) initial += bullet->getFourMomentum();
    if (target) initial += target->getFourMomentum();

    GXInuclElementaryParticle<T> const* pbullet = dynamic_cast<GXInuclElementaryParticle<T> const*>(bullet);
    GXInuclElementaryParticle<T> const* ptarget = dynamic_cast<GXInuclElementaryParticle<T> const*>(target);
    G4InuclNuclei const* nbullet = dynamic_cast<G4InuclNuclei const*>(bullet);
    G4InuclNuclei const* ntarget = dynamic_cast<G4InuclNuclei const*>(target);

    if (verboseLevel > 3) {
      std::cerr <<"======= GXBalanceCheck:\n";
      std::cerr<<"\tbullet: "<< *pbullet <<"\n";
      std::cerr<<"\ttarget: "<< *ptarget <<"\n";
      const std::vector<GXInuclElementaryParticle<T>>& outparts = output.getOutgoingParticles();
      std::cerr<<"\toutput: "<< outparts.size() <<" particles:\n";
      for(int i=0; i < outparts.size(); ++i) {
	LorentzVector<T> lorvec = outparts[i].getFourMomentum();
	std::cerr<<"  i="<< i <<" part[i]=["<< lorvec.px() <<"; "<< lorvec.py() <<"; "<< lorvec.pz() <<"; "<< lorvec.E() <<")\n";
      }
      std::cerr<<"=========\n";
    }

    // Baryon number, charge and strangeness must be computed "by hand"
    initialCharge = 0;
    if (pbullet) initialCharge += pbullet->getCharge(); // assuming homogeneous bullets and targets
    if (ptarget) initialCharge += ptarget->getCharge();

    initialBaryon =
      ((pbullet ? pbullet->baryon() : nbullet ? nbullet->getA() : 0) +
       (ptarget ? ptarget->baryon() : ntarget ? ntarget->getA() : 0) );

    // NOTE:  Currently we ignore possibility of hypernucleus target
    initialStrange = 0;
    if (pbullet) initialStrange += pbullet->getStrangeness();
    if (ptarget) initialStrange += ptarget->getStrangeness();

    // Final state totals are computed for us
    final = output.getTotalOutputMomentum();
    finalBaryon = output.getTotalBaryonNumber();
    finalCharge = output.getTotalCharge();
    finalStrange = output.getTotalStrangeness();

    // Report results
    if (verboseLevel > 2) {
      G4cout << " initial px " << initial.x() << " py " << initial.y()
	     << " pz " << initial.z() << " E " << initial.t()
	     << " baryon " << initialBaryon << " charge " << initialCharge
	     << " strange " << initialStrange << G4endl
	     << "   final px " << final.x() << " py " << final.y()
	     << " pz " << final.z() << " E " << final.t()
	     << " baryon " << finalBaryon << " charge " << finalCharge
	     << " strange " << finalStrange << G4endl;
    }
  }
  
  template <typename T>
  void GXCascadeBalanceChecker<T>::
  check(GXInuclParticle<T> const* bullet, GXInuclParticle<T> const* target,
	std::vector<GXInuclElementaryParticle<T>> const& particles)
  {
    GXCollisionOutput<T> temp;
    temp.addOutgoingParticles(particles);
    check(bullet, target, temp);
  }


  template <typename T>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE
  Mask_v<T> GXCascadeBalanceChecker<T>::
  energyOkay() const
  {
    Mask_v<T> relokay = (math::Abs(relativeE()) < relativeLimit);
    Mask_v<T> absokay = (math::Abs(deltaE()) < absoluteLimit);

    if (verboseLevel > 1 && !MaskFull(relokay | absokay)) {    // TODO: fix !relokay || !absokay
      std::cerr << fName << ": Energy conservation: relative " << relativeE()
		<< (!MaskFull(relokay) ? " conserved" : " VIOLATED")
		<< " absolute " << deltaE()
		<< (!MaskFull(absokay )? " conserved" : " VIOLATED") <<"\n";
    } else if (verboseLevel > 1) {
      std::cerr << fName << ": Energy conservation: relative " << relativeE()
		<< " conserved absolute " << deltaE() << " conserved\n";
    }

    return (relokay & absokay);
  }

  template <typename T>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE
  Mask_v<T> GXCascadeBalanceChecker<T>::
  ekinOkay() const
  {

  }

  template <typename T>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE
  Mask_v<T> GXCascadeBalanceChecker<T>::
  momentumOkay() const
  {
    Mask_v<T> relokay = (math::Abs(relativeP()) < relativeLimit);
    Mask_v<T> absokay = (math::Abs(deltaP()) < absoluteLimit);

    if (verboseLevel > 1 && !MaskFull(relokay | absokay)) {  // TODO: fix !relokay || !absokay
      std::cerr << fName << ": Momentum conservation: relative " << relativeP()
		<< (!MaskFull(relokay) ? " conserved" : " VIOLATED")
		<< " absolute " << deltaP()
		<< (!MaskFull(absokay)? " conserved" : " VIOLATED") <<"\n";
    } else if (verboseLevel > 1) {
      std::cerr << fName << ": Momentum conservation: relative " << relativeP()
		<< " conserved absolute " << deltaP() << " conserved\n";
    }

    return (relokay & absokay);
  }

  template <typename T>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE
  Mask_v<T> GXCascadeBalanceChecker<T>::
  baryonOkay() const
  {
    Mask_v<T> bokay = (vecCore::Convert<T,Int_v>(deltaB()) == T(0));	// Must be perfect!
    if (verboseLevel > 1 && !MaskFull(bokay)) {
      std::cerr << fName << ": Baryon number conservation VIOLATED " << deltaB() <<"\n";
    }
    return bokay;
  }

  template <typename T>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE
  Mask_v<T> GXCascadeBalanceChecker<T>::
  chargeOkay() const
  {
    Mask_v<T> qokay = (vecCore::Convert<T,Int_v>(deltaQ()) == T(0));	// Must be perfect!

    if (verboseLevel > 1 && !MaskFull(qokay)) {
      std::cerr << fName << ": Charge conservation VIOLATED " << deltaQ() << "\n";
    }

    return qokay;
  }

  template <typename T>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE
  Mask_v<T> GXCascadeBalanceChecker<T>::
  strangeOkay() const
  {
    Mask_v<T> sokay = (vecCore::Convert<T,Int_v>(deltaS()) == T(0));	// Must be perfect!

    if (verboseLevel && !MaskFull(sokay)) {
      std::cerr << fName << ": Strangeness conservation VIOLATED " << deltaS() <<"\n";
    }

    return sokay;
  }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXCascadeBalanceChecker_hh
