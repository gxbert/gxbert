//
// @File: GXElementaryParticleCollider.h
//
// 20180622  Guilherme Lima -- Created, based on M.Kelsey's G4ElementaryParticleCollider

#ifndef GXELEMENTARY_PARTICLE_COLLIDER_HH
#define GXELEMENTARY_PARTICLE_COLLIDER_HH

//#include "G4CascadeColliderBase.hh"
//#include "G4CascadeInterpolator.hh"
#include "GXCascadeParameters.hh"
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
#include "GXCascadeFinalStateGenerator.hh"

#include "G4CascadeChannelTables.hh"
#include "G4CascadeChannel.hh"

#include <vector>
#include <iostream>
using std::cerr;

namespace gxbert {

inline namespace GXBERT_IMPL_NAMESPACE {

// class G4ElementaryParticleCollider : public G4CascadeColliderBase {

template <typename T>
class GXElementaryParticleCollider {

  // constexpr static size_t sizeT = vecCore::VectorSize<T>();
  // using Int_v  = typename vecCore::backend::VcSimdArray<sizeT>;
  using Int_v = vecCore::Index_v<T>;
  using Bool_v = typename vecCore::Mask_v<T>;

private:
  static int verboseLevel;

  //.. Nuclear environment (to do pion-nucleon absorption)
  Index_v<T> nucleusA;
  Index_v<T> nucleusZ;

  //.. Utility class to generate final-state kinematics
  GXCascadeFinalStateGenerator<T> fsGenerator;

  //.. Internal buffers for lists of secondaries
  std::vector<GXInuclElementaryParticle<T>> particles;
  std::vector<LorentzVector<T>> scm_momentums;
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

  // Nucleus to use for recoil
  // TODO: allow distinct values here
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setNucleusState(int a, int z)
  {
    nucleusA.set( Index_v<T>(a) );
    nucleusZ.set( Index_v<T>(z) );
  }

  // Nuclei to use for recoil
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setNucleusState(Index_v<T> const& a, Index_v<T> const& z)
  {
    nucleusA.set( a );
    nucleusZ.set( z );
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void setVerboseLevel(int level) { verboseLevel = level; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void printArrays(std::vector<int> const& ikinds, Int_v mult) const;

  // Decide whether to use GXElementaryParticleCollider<T> or not
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  virtual bool useEPCollider(GXInuclParticle<T> const* bullet,
			     GXInuclParticle<T> const* target) const
  {
    return (dynamic_cast<GXInuclElementaryParticle<T> const*>(bullet) &&
	    dynamic_cast<GXInuclElementaryParticle<T> const*>(target));
  }

  void collide(GXInuclParticle<T> const* bullet,
	       GXInuclParticle<T> const* target,
	       GXCollisionOutput<T>& output);

  //private:

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Int_v generateMultiplicity(Int_v const& is, T const& ekin) const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void generateOutgoingPartTypes(Int_v is, Int_v mult, T& ekin);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void generateSCMfinalState(T& ekin, T& etot_scm,
			     GXInuclElementaryParticle<T> const* particle1,
			     GXInuclElementaryParticle<T> const* particle2);

  // Pion (or photon) absorption on a dibaryon
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void generateSCMpionAbsorption(T& etot_scm,
				 GXInuclElementaryParticle<T> const* particle1,
				 GXInuclElementaryParticle<T> const* particle2);

  // Muon absorption on a dibaryon (with outgoing neutrino)
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void generateSCMmuonAbsorption(T& etot_scm,
				 GXInuclElementaryParticle<T> const* particle1,
				 GXInuclElementaryParticle<T> const* particle2);

  /// Pion absorption on a single nucleon (charge exchange)
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void generateSCMpionNAbsorption(T& etot_scm,
				  GXInuclElementaryParticle<T> const* part1,
				  GXInuclElementaryParticle<T> const* part2,
				  vecCore::Mask<T> const& doit);

  // VECCORE_ATT_HOST_DEVICE
  // VECCORE_FORCE_INLINE
  // void rebasketizeByMultiplicity(size_t size, Index_v<T>* mult, GXInuclElementaryParticle<T>* bullets, GXInuclElementaryParticle<T>* targets) const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v pionNucleonAbsorption(T const& ekin) const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
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
      masses[i] = GXInuclElementaryParticle<T>::getParticleMass(particle_kinds[i]);
      masses2[i] = masses[i] * masses[i];
      if(verboseLevel > 3) {
	std::cerr<<"   i="<< i <<" masses: "<< masses[i] <<", masses2[i]: "<< masses2[i] <<"\n";
      }
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

/*
template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
void GXElementaryParticleCollider<T>::
rebasketizeByMultiplicity(size_t size, Index_v<T>* mult, GXInuclElementaryParticle<T>* bullets, GXInuclElementaryParticle<T>* targets) const
{
  Index_v<T> counters[10] = {0};

  // count multiplicities and build mapping
  const size_t vsize = vecCore::VectorSize<T>();
  for(size_t i = 0; i < size; i += vsize) {
    //std::cerr<<"rebask: i="<< i <<", mult="<< mult <<", bullets="<< bullets <<", targets="<< targets <<"\n";
    counters[Get(mult,i)]++;
    std::cerr<<"Loop: mult["<< i <<"]="<< mult[i] <<", counters: ";
    for(size_t j = 2; j <= 9; ++j) {
      vecCore::MaskedAssign(counters[j], mult[i] == Index_v<T>(j), counters[j] + Index_v<T>(1));
    }
  }
  std::cerr<<" counters: <";
  for (size_t j = 2; j <= 9; ++j) {
    std::cerr<<"["<< j <<':'<< counters[j]<<"], ";
  }
  std::cerr<<"\n";

  // swap tracks as appropriate
}
*/

template <typename T>
vecCore::Index_v<T> GXElementaryParticleCollider<T>::
generateMultiplicity(vecCore::Index_v<T> const& hadPairs, T const& ekin) const
{
  Int_v mult(0);
  static bool first = true;

  const G4CascadeChannel* xsecTable;
  size_t vsize = vecCore::VectorSize<T>();
  int lastPair = -999;
  //std::cerr<<" genMult(): vsize="<< vsize <<", lastPair="<< lastPair <<", hadPairs="<< hadPairs <<" and ekin="<< ekin <<"\n";
  for (size_t i =0; i < vsize; ++i) {
    double laneEkin = Get(ekin, i);
    int lanePair    = Get(hadPairs, i);
    //std::cerr<<" genMult(): i="<< i <<", laneEkin="<< laneEkin <<", lanePair="<< lanePair <<" lastPair="<< lastPair <<"\n";
    if (lanePair != lastPair) {
      xsecTable = G4CascadeChannelTables::GetTable(lanePair);
      lastPair = lanePair;
    }

    //std::cerr<<" genMult(): xsecTable="<< xsecTable <<", lastPair="<< lastPair <<"\n";
    if(first) {
      //xsecTable->printTable(std::cerr);
      first = false;
    }
    if (xsecTable) {
      int temp = xsecTable->getMultiplicity(laneEkin);
      Set(mult, i, temp);
    }
    else {
      std::cerr << " GXElementaryParticleCollider: Unknown interaction channel "
	<< lanePair << " - multiplicity not generated.\n";
    }
  }

  if(verboseLevel > 3) {
    std::cerr << " GXElementaryParticleCollider::generateMultiplicity: "
      << " multiplicity = " << mult <<"\n";
  }

  //return mult;
  ////// temporary!!! return fixed multiplicity
  return Index_v<T>(4);
}

/* // GL Note: I tried to push vectorization into xsecTable->GetMultiplicity(ekin), but got "virtual template methods not allowed" problems
template <typename T>
Index_v<T> GXElementaryParticleCollider<T>::generateMultiplicity(Int_v const& hadPairs, T const& ekin) const
{
  Int_v zero(0);
  Int_v mul = zero;

  int i = 0;
  size_t vsize = vecCore::VectorSize<T>();
  while ( i < vsize || vecCore::MaskFull(mul > 0) ) {
    int nextLane = Find(mul, 0);         // find lane where multiplicity has not yet been drawn
    int nextPair = Get(hadPairs, nextLane);  // get hadron pair in that lane

    const GXCascadeChannel* xsecTable = GXCascadeChannelTables::GetTable(nextPair);
    if (xsecTable) {
      const Int_v temp = xsecTable->getMultiplicity(ekin);
      vecCore::MaskedAssign( mul, hadPairs == nextPair, temp); // update all lanes as appropriate
    }
    else {
      std::cerr << " GXElementaryParticleCollider: Unknown interaction channel "
	<< nextPair << " - multiplicity not generated.\n";
    }
    i++;
  }

  if(verboseLevel > 3){
    std::cerr << " GXElementaryParticleCollider::generateMultiplicity: "
    << " multiplicity = " << mul <<"\n";
  }

  return mul;
}
*/

template <typename T>
void GXElementaryParticleCollider<T>::
generateOutgoingPartTypes(Int_v hadPairs, Int_v mult, T& ekin)
{
  assert(isHomogeneous(hadPairs) && "Non-homogeneous input initial state.");
  if (!isHomogeneous(hadPairs)) {
    std::cerr<< ">>> generateOutgoingPartTypes(): ERROR: Non-homogeneous input initial state: "<< hadPairs <<"\n";
  }

  assert(isHomogeneous(mult) && "Non-homogeneous input multiplicity.");
  if (!isHomogeneous(mult)) {
    std::cerr<< ">>> generateOutgoingPartTypes(): ERROR: Non-homogeneous input initial state: "<< hadPairs <<"\n";
  }

  particle_kinds.clear();	// Initialize buffer for generation

  //=== original code, for reference -- not hard to vectorize if given homogeneous input
//   int initState = Get(hdPairs, 0);
//   const G4CascadeChannel* xsecTable = G4CascadeChannelTables::GetTable(initState);
//   if (xsecTable)
//     xsecTable->getOutgoingParticleTypes(particle_kinds, mult, ekin);
//   else {
//     std::cerr << " GXElementaryParticleCollider: Unknown interaction channel "
//       << initState << " - outgoing kinds not generated.\n";
//   }

//   return;
// }
  //=== end of original code

  const size_t cmultiplicity = Get(mult, 0);
  particle_kinds.resize(cmultiplicity);
  std::vector<int> tempKinds;  // buffer for scalar calls to xsecTable->getOutgoingParticleTypes()

  const G4CascadeChannel* xsecTable;
  size_t vsize = vecCore::VectorSize<T>();
  int lastPair = -999;
  for (size_t i =0; i < vsize; ++i) {
    double laneEkin = Get(ekin, i);
    int lanePair    = Get(hadPairs, i);
    size_t nkinds = particle_kinds.size();
    if (lanePair != lastPair) {
      xsecTable = G4CascadeChannelTables::GetTable(lanePair);
      lastPair = lanePair;
    }
    //std::cerr<<" i="<< i <<", laneEkin="<< laneEkin <<", lanePair="<< lanePair <<", lastPair="<< lastPair <<", xsecTable="<<  xsecTable <<"\n";
    if (xsecTable) {
      xsecTable->getOutgoingParticleTypes(tempKinds, nkinds, laneEkin);
      // std::cerr<<" xsecTable: tempKinds["<< tempKinds.size() <<"] = [";
      // for (size_t ii=0; ii < cmultiplicity; ++ii) {
      //   std::cerr<< tempKinds[ii] <<", ";
      // }
      // std::cerr<<"]\n";
    }
    else {
      std::cerr << " GXElementaryParticleCollider: Unknown interaction channel "
	<< lanePair << " - outgoing kinds not generated.\n";
    }

    for (size_t isec = 0; isec < cmultiplicity; ++isec) {
      Set( particle_kinds[isec], i, tempKinds[isec] );
    }
  }

  if (verboseLevel > 3) {
    std::cerr<<" returning from generateOutPartTypes: "; this->printArrays(tempKinds, mult);
  }

  return;
}

template <typename T>
void GXElementaryParticleCollider<T>::printArrays(std::vector<int> const& ikinds, Int_v mult) const
{
  const size_t cmultiplicity = Get(mult, 0);
  const size_t ksize = ikinds.size();
  std::cerr<<" mult="<< mult
	   <<", ikinds[size="<< ksize <<"]: [";
  for (size_t i = 0; i < ksize; ++i) {
    std::cerr << ikinds[i] <<" ";
  }
  std::cerr<<"] and part_kinds: [";
  for (size_t i = 0; i < cmultiplicity; ++i) {
    std::cerr << particle_kinds[i] <<" ";
  }
  std::cerr <<"]\n";
}

template <typename T>
void GXElementaryParticleCollider<T>::generateSCMfinalState(T& ekin, T& etot_scm,
							    GXInuclElementaryParticle<T> const* particle1,
							    GXInuclElementaryParticle<T> const* particle2)
{
  constexpr int itry_max = 10;

  if (verboseLevel > 2) {
    std::cerr << ">>> GXElementaryParticleCollider::generateSCMfinalState()...\n";
  }

  fsGenerator.SetVerboseLevel(verboseLevel);

  Int_v type1 = particle1->type();
  Int_v type2 = particle2->type();

  Int_v is = type1 * type2;
  if (verboseLevel > 3) std::cerr << " is " << is <<"\n";

  Int_v multiplicity = 0;
  bool goodStates = false;

  // Initialize buffers for this event
  particles.clear();
  particle_kinds.clear();

  int itry = 0;
  // Generate list of final-state particles
  multiplicity = generateMultiplicity(is, ekin);

  //.. rebasketize cf. multiplicity
  //... not needed for now, since multiplicity was forced to be constant in generateMultiplicity()

  // no need to be inside the loop (perf improvement)
  generateOutgoingPartTypes(is, multiplicity, ekin);
  if (particle_kinds.empty()) {
    if (verboseLevel > 3) {
      std::cerr << " generateOutgoingPartTypes failed mult " << multiplicity <<"\n";
    }
    //continue;
  }

  fillOutgoingMasses();	// Fill mass buffer from particle types

  fsGenerator.Configure(particle1, particle2, particle_kinds);

  while (!goodStates && itry++ < itry_max) {
    // else {  // debugging printout only!
    //   cerr<<" generateOutgoingPartTypes(): partKinds=[";
    //   typename std::vector<Int_v>::const_iterator iter = particle_kinds.begin();
    //   for( ; iter != particle_kinds.end(); ++iter ) {
    // 	std::cerr<< *iter <<' ';
    //   }
    //   std::cerr<<"]\n";
    // }

    goodStates = fsGenerator.Generate(etot_scm, masses, scm_momentums);
  }	// while (generate) 

  if (itry >= itry_max) {		// Unable to generate valid final state
    if (verboseLevel > 2)
      std::cerr << " generateSCMfinalState failed " << itry << " attempts\n";
    return;
  }

  // Store generated momenta into outgoing particles

  size_t maxmult = vecCore::ReduceMax(multiplicity);
  particles.resize(maxmult);		// Preallocate buffer
  //std::cerr<<" ***** spot 1 - particles.size="<< particles.size() <<", maxmult="<< maxmult <<"\n";
  for (size_t i=0; i<maxmult; i++) {
    //std::cerr<<" i="<< i <<", scm_mom[i]="<< scm_momentums[i] <<", kinds[i]="<< particle_kinds[i] <<"\n";
    particles[i].fill(scm_momentums[i], particle_kinds[i], gxbert::EPCollider);
  }

  if (verboseLevel > 3) {
    std::cerr << " <<< GXElementaryParticleCollider::generateSCMfinalState.\n";
  }

  return;	// Particles buffer filled
}

// Pion (or photon) absorption on a dibaryon
template <typename T>
void GXElementaryParticleCollider<T>::
generateSCMpionAbsorption(T& etot_scm,
			  GXInuclElementaryParticle<T> const* particle1,
			  GXInuclElementaryParticle<T> const* particle2)
{

}

// Muon absorption on a dibaryon (with outgoing neutrino)
template <typename T>
void GXElementaryParticleCollider<T>::
generateSCMmuonAbsorption(T& etot_scm,
			  GXInuclElementaryParticle<T> const* particle1,
			  GXInuclElementaryParticle<T> const* particle2)
{

}

template <typename T>
void GXElementaryParticleCollider<T>::
generateSCMpionNAbsorption(T& etot_scm,
			   GXInuclElementaryParticle<T> const* part1,
			   GXInuclElementaryParticle<T> const* part2,
			   vecCore::Mask<T> const& doit)
{
  if (verboseLevel > 3)
    cerr << " >>> GXElementaryParticleCollider::generateSCMpionNAbsorption\n";

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
  Int_v ntype(type1);
  MaskedAssign(ntype, part2->nucleon(), type2);
  Int_v outType = Int_v(3) - ntype;
  vecCore::MaskedAssign(outType, !doit, Int_v(-1));
  particle_kinds.push_back(outType);

  fillOutgoingMasses();

  // Get mass of residual nucleus (2-ntype = 1 for proton, 0 for neutron)
  T mRecoil(0.);
  for (size_t i = 0; i < VectorSize<T>(); ++i) {
    double newA = Get(nucleusA,i) - 1;
    double newZ = Get(nucleusZ,i) - (2 - Get(ntype,i));
    Set(mRecoil, i, G4InuclNuclei::getNucleiMass(newA, newZ));
  }
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
  LorentzVector<T> mom1 = GXInuclSpecialFunctions::generateWithRandomAnglesNew<T>(pmod, masses[0]);
  //LorentzVector<T> mom1 = generateWithRandomAnglesNew<T>(pmod, masses[0]);  // compilation error: no template named 'gener...glesNew'

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
  particles[0].fill(mom1, particle_kinds[0], EPCollider);
}

template <typename T>
void GXElementaryParticleCollider<T>::
collide(GXInuclParticle<T> const* bullet, GXInuclParticle<T> const* target, GXCollisionOutput<T>& output)
{
  if (verboseLevel > 1) cerr << " >>> GXElementaryParticleCollider::collide\n";

  // Sanity check
  Bool_v done(!useEPCollider(bullet,target));
  if ( !vecCore::MaskEmpty(done) ) {
    cerr << " ElementaryParticleCollider -> can collide only particle with particle: done="<< done <<"\n";
    return;
  }

#ifdef G4CASCADE_DEBUG_SAMPLER
  static bool doPrintTables = true;	// Once and only once per job
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

  if (!particle1 || !particle2) {    // Redundant with useEPCollider() ?
    cerr << " ElementaryParticleCollider -> can only collide hadrons!  done="<< done <<"\n";
    return;
  }

  // Keep track of lanes for which no cascading is needed
  const size_t vsize = vecCore::VectorSize<T>();
  done = done | Bool_v( particle1->isNeutrino() || particle2->isNeutrino() );

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
    Bool_v doPionNAbsorption = pionNucleonAbsorption(ekin);
    //cerr<<" doPionNAbsorption="<< doPionNAbsorption <<"\n";

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
    // for now, use the existing scalar function -- ok since particle1,2 are homogeneous!
    int ptype1 = Get(particle1->type(), 0);
    int ptype2 = Get(particle2->type(), 0);
    if (!G4NucleiModel::useQuasiDeuteron(ptype1, ptype2) &&
	!G4NucleiModel::useQuasiDeuteron(ptype2, ptype1)) {
      cerr << " ElementaryParticleCollider -> can only collide pi,mu,gamma with dibaryons\n";
      return;
    }

    if (MaskFull(particle1->isMuon()) || MaskFull(particle2->isMuon())) {
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

  //   cerr<< " <<< Non-conservation in GXElementaryParticleCollider\n";
  // }

  // in vector mode, sorting is painful!
  //std::sort(particles.begin(), particles.end(), GXParticleLargerEkin<double>());
  output.addOutgoingParticles(particles);
  if (verboseLevel > 1) {
    std::cerr <<" GXEPCollider::collide(): returning... particle[0]="<< particles[0] <<" and part[1]="<< particles[1] << G4endl;
  }
}

template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
vecCore::Mask_v<T> GXElementaryParticleCollider<T>::
pionNucleonAbsorption(T const& ekin) const
{
  if (verboseLevel > 3)
    cerr << " >>> GXElementaryParticleCollider::pionNucleonAbsorption ?"
	 << " ekin " << ekin << " is " << interCase.hadrons() <<"\n";

  // Absorption occurs with specified probability
  const T absProb(GXCascadeParameters::piNAbsorption());

  // Absorption occurs only for pi- p -> n, or pi+ n -> p
  // Restrict to "very slow" pions, to allow for some normal scattering
  return ((interCase.hadrons() == pim*pro || interCase.hadrons() == pip*neu)
	  && (ekin < 0.05)		// 50 MeV kinetic energy or less
	  && (inuclRndm() < absProb)
	  );
}

/*
// saves final state particles into a SIMD array similar to the input format
template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
void GXElementaryParticleCollider<T>::
SaveOutgoingParticles(G4CollisionOutput& output) const
{
  std::vector<GXInuclElementaryParticle<T>>::const_iterator ipartSIMD = particles.begin();
  for( ; ipartSIMD != particles.end(); ++iparticles) {

  }
}
*/

// initialize verboseLevel
template <typename T>
int GXElementaryParticleCollider<T>::verboseLevel = 0;

}
}
#endif	// GXELEMENTARY_PARTICLE_COLLIDER_HH
