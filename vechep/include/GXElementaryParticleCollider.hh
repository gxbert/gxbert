//
// @File: GXElementaryParticleCollider.h
//
// 20180622  Guilherme Lima -- Created, based on M.Kelsey's G4ElementaryParticleCollider

#ifndef GXELEMENTARY_PARTICLE_COLLIDER_HH
#define GXELEMENTARY_PARTICLE_COLLIDER_HH

#include "G4CascadeColliderBase.hh"
#include "G4CascadeFinalStateGenerator.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeParameters.hh"
#include "G4NucleiModel.hh"
#include "LorentzVector.hh"

#include "GXInuclParticle.hh"
#include "GXInuclElementaryParticle.hh"
#include "GXInteractionCase.hh"
#include "LorentzVector.hh"
#include "GXLorentzConvertor.hh"
#include "GXParticleLargerEkin.hh"
#include "GXCollisionOutput.hh"
#include "GXInuclSpecialFunctions.hh"

#include <vector>
#include <iostream>
using std::cerr;

namespace gxbert {

inline namespace GXBERT_IMPL_NAMESPACE {

// class G4ElementaryParticleCollider : public G4CascadeColliderBase {

template <typename T>
class GXElementaryParticleCollider {

  using Int_v  = vecCore::Index_v<T>;
  using Bool_v = vecCore::Mask_v<T>;

private:
  static int verboseLevel;

  //.. Nuclear environment (to do pion-nucleon absorption)
  Index_v<T> nucleusA;
  Index_v<T> nucleusZ;

  //.. Utility class to generate final-state kinematics
  G4CascadeFinalStateGenerator fsGenerator;

  //.. Internal buffers for lists of secondaries
  std::vector<GXInuclElementaryParticle<T>> particles;
  //std::vector<LorentzVector<T>> scm_momenta;
  std::vector<T> masses;
  std::vector<T> masses2;
  std::vector<Int_v> particle_kinds;

  //.. used to be on a base class
  GXInteractionCase<T> interCase;

public:
  GXElementaryParticleCollider()
    : nucleusA(0)
    , nucleusZ(0)
  { }

  virtual ~GXElementaryParticleCollider()
  { }

  void collide(GXInuclParticle<T> const* bullet, GXInuclParticle<T> const* target, GXCollisionOutput& output)
  {
    if (verboseLevel > 1)
      cerr << " >>> G4ElementaryParticleCollider::collide\n";

    // Sanity check
    Bool_v done = !useEPCollider(bullet,target);
    if ( !vecCore::MaskEmpty(done) ) {
      cerr << " ElementaryParticleCollider -> can collide only particle with particle: done="<< done <<"\n";
      return;
    }

#ifdef G4CASCADE_DEBUG_SAMPLER
    static G4bool doPrintTables = true;	// Once and only once per job
    if (doPrintTables) {
      G4CascadeChannelTables::Print();
      doPrintTables = false;
    }
#endif

    interCase.set(bullet, target);	// To identify kind of collision

    GXInuclElementaryParticle<T> const* particle1 = dynamic_cast<GXInuclElementaryParticle<T> const*>(bullet);
    GXInuclElementaryParticle<T> const* particle2 = dynamic_cast<GXInuclElementaryParticle<T> const*>(target);

    if (verboseLevel > 1) {
      cerr <<" GXEPCollider::collide() - input particles:"<< *particle1 << "\n and: " << *particle2 <<"\n";
    }

    done = done | (!particle1) | (!particle2);
    //if (!particle1 || !particle2) {    // Redundant with useEPCollider()
    if ( !vecCore::MaskEmpty(done) ) {    // Redundant with useEPCollider()
      cerr << " ElementaryParticleCollider -> can only collide hadrons!  done="<< done <<"\n";
      return;
    }

    // Keep track of lanes for which no cascading is needed
    const size_t vsize = vecCore::VectorSize<T>();
    done = done | ( particle1->isNeutrino() || particle2->isNeutrino() );

    // Check if input is homogeneous -- redundant!?
    Index_v<T> hadcase = interCase.hadrons();
    const size_t case0 = Get(hadcase, 0);
    bool allsame(true);
    for (size_t i = 1; i < vsize; ++i) {
      allsame &= (case0 == Get(hadcase, i));
    }
    assert(allsame && "GXEPCollider::collide() called for inhomogeneous interCases");

    // Check for available interaction table, if not quasi-deuteron special cases
    {
      const Bool_v cantCollide = (!particle1->quasi_deutron() && !particle2->quasi_deutron() && !interCase.hasValidTable());
      done = done | cantCollide;
      if (!vecCore::MaskEmpty(cantCollide)) {
	cerr << " ElementaryParticleCollider -> cantCollide="<< cantCollide <<" for "
	     << *particle1 << " and "<< *particle2 <<" - done="<< done <<"\n";
      }
    }

    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    GXLorentzConvertor<T> convertToSCM;    // Utility to handle frame manipulation
    convertToSCM.setVerbose(verboseLevel);

    // TODO: need to swap appropriately at this point!!!
    if (vecCore::MaskFull(particle2->nucleon() | particle2->quasi_deutron()) ) {
      convertToSCM.setBullet(*particle1);
      convertToSCM.setTarget(*particle2);
    } else {
      cerr<<" ***** GXEPCollider::collide(): SWAP WAS NEEDED!!  CHECK RESULTS!!\n";
      convertToSCM.setBullet(*particle2);
      convertToSCM.setTarget(*particle1);
    }
    convertToSCM.toTheCenterOfMass();

    T etot_scm = convertToSCM.getTotalSCMEnergy();

    // Generate any particle collision with nucleon
    if (MaskFull(particle1->nucleon()) || MaskFull(particle2->nucleon())) {
      T ekin = convertToSCM.getKinEnergyInTheTRS();

      // SPECIAL: very low energy pions may be absorbed by a nucleon
      Bool_v doPionNAbsorption = vecCore::MaskEmpty(pionNucleonAbsorption(ekin));
      cerr<<" doPionNAbsorption="<< doPionNAbsorption <<"\n";

      // use if() here to save some time, since will be mostly false
      if ( !vecCore::MaskEmpty(doPionNAbsorption) ) {
	generateSCMpionNAbsorption(etot_scm, particle1, particle2, doPionNAbsorption);
	done = done | doPionNAbsorption;
      }

      if ( !vecCore::MaskFull(doPionNAbsorption) ) {
	generateSCMfinalState(ekin, etot_scm, particle1, particle2);
      }
    }

    // Generate pion or photon collision with quasi-deuteron
    if (MaskFull(particle1->quasi_deutron()) || MaskFull(particle2->quasi_deutron())) {
      if (!G4NucleiModel::useQuasiDeuteron(particle1->type(), particle2->type()) &&
	  !G4NucleiModel::useQuasiDeuteron(particle2->type(), particle1->type())) {
	cerr << " ElementaryParticleCollider -> can only collide pi,mu,gamma with dibaryons\n";
	return;
      }

      if (particle1->isMuon() || particle2->isMuon()) {
	generateSCMmuonAbsorption(etot_scm, particle1, particle2);
      } else {		// Currently, pion absoprtion also handles gammas
	generateSCMpionAbsorption(etot_scm, particle1, particle2);
      }
    }

    if (particles.empty()) {	// No final state possible, pass bullet through
      if (verboseLevel) {
	cerr << " ElementaryParticleCollider -> failed to collide "
	     << particle1->getMomModule() << " GeV/c "
	     << *particle1 << " with "
	     << *particle2 <<"\n";
      }
      return;
    }

    // Convert final state back to lab frame
    // TODO: vectorize this loop?
    LorentzVector<T> mom;      	// Buffer to avoid memory churn
    typename std::vector<GXInuclElementaryParticle<T>>::iterator ipart = particles.begin();
    for( ; ipart != particles.end(); ipart++) {	
      mom = convertToSCM.backToTheLab(ipart->getFourMomentum());
      ipart->setMomentum(mom); 
    }
  
    // // Check conservation in multibody final state
    // if (verboseLevel && !validateOutput(bullet, target, particles)) {
    //   cerr << " incoming particles: \n" << *particle1 << G4endl
    // 	   << *particle2 <<"\n"
    // 	   << " outgoing particles:\n";
    //   for(ipart = particles.begin(); ipart != particles.end(); ipart++)
    // 	cerr << *ipart <<"\n";

    //   cerr<< " <<< Non-conservation in G4ElementaryParticleCollider\n";
    // }
    
    std::sort(particles.begin(), particles.end(), GXParticleLargerEkin<double>());
    output.addOutgoingParticles(particles);
  }

  // TODO: allow distinct values here
  void setNucleusState(int a, int z) {	// Nucleus to use for recoil
    nucleusA.set( Index_v<T>(a) );
    nucleusZ.set( Index_v<T>(z) );
  }

  void setNucleusState(Index_v<T> const& a, Index_v<T> const& z) {  // Nuclei to use for recoil
    nucleusA.set( a );
    nucleusZ.set( z );
  }

  void setVerboseLevel(int level) { verboseLevel = level; }

  // Decide whether to use G4ElementaryParticleCollider or not
  virtual vecCore::Mask_v<T> useEPCollider(GXInuclParticle<T> const* bullet,
					   GXInuclParticle<T> const* target) const
  {
    return (dynamic_cast<GXInuclElementaryParticle<T> const*>(bullet) &&
	    dynamic_cast<GXInuclElementaryParticle<T> const*>(target));
  }


private:

  Int_v generateMultiplicity(Int_v const& is, T const& ekin) const;

  void generateOutgoingPartTypes(int is, int mult, T& ekin);

  void generateSCMfinalState(T& ekin, T& etot_scm,
			     GXInuclElementaryParticle<T> const* particle1,
			     GXInuclElementaryParticle<T> const* particle2);

  // Pion (or photon) absorption on a dibaryon
  void generateSCMpionAbsorption(T& etot_scm,
				 GXInuclElementaryParticle<T> const* particle1,
				 GXInuclElementaryParticle<T> const* particle2);

  // Muon absorption on a dibaryon (with outgoing neutrino)
  void generateSCMmuonAbsorption(T& etot_scm,
				 GXInuclElementaryParticle<T> const* particle1,
				 GXInuclElementaryParticle<T> const* particle2);

  /// Pion absorption on a single nucleon (charge exchange)
  void generateSCMpionNAbsorption(T& etot_scm,
				  GXInuclElementaryParticle<T> const* part1,
				  GXInuclElementaryParticle<T> const* part2,
				  vecCore::Mask<T> const& doit);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v pionNucleonAbsorption(T const& ekin) const
  {
    if (verboseLevel > 3)
      cerr << " >>> GXElementaryParticleCollider::pionNucleonAbsorption ?"
	   << " ekin " << ekin << " is " << interCase.hadrons() <<"\n";

    // Absorption occurs with specified probability
    const T absProb(G4CascadeParameters::piNAbsorption());

    // Absorption occurs only for pi- p -> n, or pi+ n -> p
    // Restrict to "very slow" pions, to allow for some normal scattering
    return ((interCase.hadrons() == pim*pro || interCase.hadrons() == pip*neu)
	    && (ekin < 0.05)		// 50 MeV kinetic energy or less
	    && (inuclRndm() < absProb)
	    );
  }

  bool splitQuasiDeuteron(int qdtype); 	// Fill kinds with NN components

  //.. Fill mass arrays from particle types
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void fillOutgoingMasses()
  {
    size_t mult = particle_kinds.size();
    masses.resize(mult,0.);
    masses2.resize(mult,0.);		// Allows direct [i] setting

    for (size_t i = 0; i < mult; i++) {
      // uses GXInuclElemParticle::setParticleMass()
      particles[i].setParticleMass(particle_kinds[i]);
      masses[i] = particles[0].mass();
      masses2[i] = masses[i] * masses[i];
    }
  }

private:
  // Copying of modules is forbidden
  GXElementaryParticleCollider(const GXElementaryParticleCollider&);
  GXElementaryParticleCollider& operator=(const GXElementaryParticleCollider&);

  // temporary!!!
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T inuclRndm() const {
    static T val(0.00234);

    val += 0.0001342;
    vecCore::MaskedAssign( val, val > T(1.0), val - T(1.0));

    return val;
  }
};

template <typename T>
Index_v<T> GXElementaryParticleCollider<T>::generateMultiplicity(Int_v const& hadPairs, T const& ekin) const
{/*
  Int_v zero(0);
  Int_v mul = zero;

  int i = 0;
  size_t vsize = vecCore::VectorSize<T>();
  while ( i < vsize || vecCore::MaskFull(mul > 0) ) {
    if Get(mul, i);
    int nextLane = Find(mul, 0);         // find lane where multiplicity has not yet been drawn
    int nextPair = Get(hadPairs, next);  // get hadron pair in that lane

    const G4CascadeChannel* xsecTable = G4CascadeChannelTables::GetTable(nextPair);
    if (xsecTable) {
      const Int_v temp = xsecTable->getMultiplicity(ekin);
      vecCore::MaskedAssign( mul, hadPairs == nextPair, temp); // update all lanes as appropriate
    }
    else {
      std::cerr << " G4ElementaryParticleCollider: Unknown interaction channel "
	<< nextPair << " - multiplicity not generated.\n";
    }
    i++;
  }

  if(verboseLevel > 3){
    G4cout << " G4ElementaryParticleCollider::generateMultiplicity: "  
           << " multiplicity = " << mul << G4endl; 
  }

  return mul;
 */
}

template <typename T>
void GXElementaryParticleCollider<T>::generateSCMfinalState(T& ekin, T& etot_scm,
							    GXInuclElementaryParticle<T> const* particle1,
							    GXInuclElementaryParticle<T> const* particle2)
{
  /*
  constexpr int itry_max = 10;

  if (verboseLevel > 2) {
    std::cerr << ">>> G4ElementaryParticleCollider::generateSCMfinalState\n";
  }

  fsGenerator.SetVerboseLevel(verboseLevel);

  Int_v type1 = particle1->type();
  Int_v type2 = particle2->type();

  Int_v is = type1 * type2;
  if (verboseLevel > 3) std::cerr << " is " << is <<"\n";

  Int_v multiplicity = 0;
  Bool_v generate    = true;

  // Initialize buffers for this event
  particles.clear();
  particle_kinds.clear();

  int itry = 0;
  while (generate && itry++ < itry_max) {

    // Generate list of final-state particles
    multiplicity = generateMultiplicity(is, ekin);

    generateOutgoingPartTypes(is, multiplicity, ekin);
    if (particle_kinds.empty()) {
      if (verboseLevel > 3) {
	G4cout << " generateOutgoingPartTypes failed mult " << multiplicity
	       << G4endl;
      }
      continue;
    }

    fillOutgoingMasses();	// Fill mass buffer from particle types

    // Attempt to produce final state kinematics
    fsGenerator.Configure(particle1, particle2, particle_kinds);
    generate = !fsGenerator.Generate(etot_scm, masses, scm_momentums);
  }	// while (generate) 

  if (itry >= itry_max) {		// Unable to generate valid final state
    if (verboseLevel > 2)
      G4cout << " generateSCMfinalState failed " << itry << " attempts"
	     << G4endl;
    return;
  }

  // Store generated momenta into outgoing particles

  particles.resize(multiplicity);		// Preallocate buffer
  for (G4int i=0; i<multiplicity; i++) {
    particles[i].fill(scm_momentums[i], particle_kinds[i],
		      G4InuclParticle::EPCollider);
  }

  if (verboseLevel > 3) {
    G4cout << " <<< G4ElementaryParticleCollider::generateSCMfinalState"
	   << G4endl;
  }
  */
  return;	// Particles buffer filled
}

// Pion (or photon) absorption on a dibaryon
template <typename T>
void GXElementaryParticleCollider<T>::generateSCMpionAbsorption(T& etot_scm,
								GXInuclElementaryParticle<T> const* particle1,
								GXInuclElementaryParticle<T> const* particle2)
{

}

// Muon absorption on a dibaryon (with outgoing neutrino)
template <typename T>
void GXElementaryParticleCollider<T>::generateSCMmuonAbsorption(T& etot_scm,
								GXInuclElementaryParticle<T> const* particle1,
								GXInuclElementaryParticle<T> const* particle2)
{

}

template <typename T>
void GXElementaryParticleCollider<T>::generateSCMpionNAbsorption(T& etot_scm,
								 GXInuclElementaryParticle<T> const* part1,
								 GXInuclElementaryParticle<T> const* part2,
								 vecCore::Mask<T> const& doit)
{
  if (verboseLevel > 3)
    cerr << " >>> G4ElementaryParticleCollider::generateSCMpionNAbsorption\n";

  particles.clear();		// Initialize buffers for this event
  particles.resize(1);

  particle_kinds.clear();

  const Int_v &type1 = part1->type();
  const Int_v &type2 = part2->type();

  // Ensure that single-nucleon absportion is valid (charge exchangeable)
  Bool_v skip = (type1*type2 != pim*pro && type1*type2 != pip*neu);
  if ( !vecCore::MaskEmpty(skip & doit) ) {
    cerr << " *** Problems?! GXElemPartCollider::generateSCMpionNAbsorption(): " << *part1 << " + " << *part2 << " -> ? - doit="<< doit
	 <<" and skip="<< skip <<"\n";
    //return;
  }

  // Get outgoing nucleon type using charge exchange
  // Proton code is 1, neutron code is 2, so 3-# swaps them
  const Int_v ntype = (part2->nucleon() ? type2 : type1);
  Int_v outType = 3 - ntype;
  vecCore::MaskedAssign(outType, !doit, Int_v(-1));
  particle_kinds.push_back(outType);

  fillOutgoingMasses();

  // Get mass of residual nucleus (2-ntype = 1 for proton, 0 for neutron)
  T mRecoil = G4InuclNuclei::getNucleiMass(nucleusA - Int_v(1), nucleusZ - (Int_v(2) - ntype));
  T mRecoil2 = mRecoil * mRecoil;

  // Recompute Ecm to include nucleus (for recoil kinematics)
  LorentzVector<T> piN4 = part1->getFourMomentum() + part2->getFourMomentum();
  T zero(0.0);
  LorentzVector<T> vsum(zero, zero, zero, mRecoil);
  vsum += piN4;

  // Two-body kinematics (nucleon against nucleus) in overall CM system
  T esq_scm = vsum.M2();
  T a = T(0.5) * (esq_scm - masses2[0] - mRecoil2);

  T pmod = vecCore::math::Sqrt((a * a - masses2[0] * mRecoil2) / esq_scm);
  //LorentzVector<T> mom1 = gxbert::GXInuclSpecialFunctions::generateWithRandomAngles<VectorBackend>(pmod, masses[0]);
  LorentzVector<T> mom1 = gxbert::GXInuclSpecialFunctions::generateWithRandomAnglesNew<T>(pmod, masses[0]);

  if (verboseLevel > 3) {
    cerr << " outgoing type " << outType << " recoiling on nuclear mass "
	 << mRecoil << "\n a " << a << " p " << pmod << " Ekin "
	 << mom1.E() - masses[0] << G4endl;
  }

  mom1.Boost(-piN4.BoostVector());	// Boost into CM of pi-N collision

  if (verboseLevel > 3) {
    cerr << " in original pi-N frame p(SCM) " << mom1.Mag() << " Ekin "
	 << mom1.E() - masses[0] <<"\n";
  }

  // Fill only the ejected nucleon
  particles[0].fill(mom1, particle_kinds[0], gxbert::EPCollider);
}

// initialize verboseLevel
template <typename T>
int GXElementaryParticleCollider<T>::verboseLevel = 0;

}
}
#endif	// GXELEMENTARY_PARTICLE_COLLIDER_HH
