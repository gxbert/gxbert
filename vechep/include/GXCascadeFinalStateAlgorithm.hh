//
// @File: GXCascadeFinalStateAlgorithm.hh
//
// 20180713 Guilherme Lima -- Created, based on the original M.Kelsey's G4CascadeFinalStateAlgorithm

#ifndef GXCascadeFinalStateAlgorithm_hh
#define GXCascadeFinalStateAlgorithm_hh 1

#include "GXVCascadeAlgorithm.hh"
#include "GXLorentzConvertor.hh"
#include "GXInuclSpecialFunctions.hh"
#include "GXTwoBodyAngularDist.hh"
#include "GXMultiBodyMomentumDist.hh"
#include "GXPowVec.hh"

#include "G4VTwoBodyAngDst.hh"
//#include "G4VMultiBodyMomDst.hh"

#include <iostream>
#include <vector>
#include <numeric>

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

class G4VMultiBodyMomentumDist;

using GXInuclSpecialFunctions::inuclRndm;

template <typename T> class GXInuclElementaryParticle;

constexpr int itry_max = 10;		// Maximum number of generation attempts

template <typename T>
class GXCascadeFinalStateAlgorithm : public GXVCascadeAlgorithm<T> {

  using Bool_v = vecCore::Mask_v<T>;
  using GXVCascadeAlgorithm<T>::GXVCascadeAlgorithm;
  using GXVCascadeAlgorithm<T>::GetName;
  using GXVCascadeAlgorithm<T>::GetVerboseLevel;

public:
  GXCascadeFinalStateAlgorithm();
  virtual ~GXCascadeFinalStateAlgorithm();

  virtual void SetVerboseLevel(int verbose);  // Pass through to factories

  // Select appropriate distributions based on interaction
  void Configure(GXInuclElementaryParticle<T>* bullet,
		 GXInuclElementaryParticle<T>* target,
		 const std::vector<Int_v>& particle_kinds);

protected:
  // Two-body generation uses angular-distribution function
  virtual void GenerateTwoBody(T initialMass,
			       std::vector<T> const& masses,
			       std::vector<LorentzVector<T>>& finalState);

  // N-body generation uses momentum-modulus distribution, computed angles
  virtual void GenerateMultiBody(T initialMass,
				 const std::vector<T>& masses,
				 std::vector<LorentzVector<T>>& finalState);

  // Compute kinematic quantities needed for distributions
  void SaveKinematics(GXInuclElementaryParticle<T>* bullet,
		      GXInuclElementaryParticle<T>* target);

  // Select generator based on initial and final state
  void ChooseGenerators(Int_v is, Int_v fs);

  // Generate momentum magnitudes and validate for use
  void FillMagnitudes(T initialMass,
		      const std::vector<T>& masses);

  vecCore::Mask_v<T> satisfyTriangle(const std::vector<T>& pmod) const;

  // Generate momentum directions into final state
  void FillDirections(T initialMass,
		      const std::vector<T>& masses,
		      std::vector<LorentzVector<T>>& finalState);

  void FillDirThreeBody(T initialMass,
			const std::vector<T>& masses,
			std::vector<LorentzVector<T>>& finalState, Bool_v mult3);

  void FillDirManyBody(T initialMass,
		       const std::vector<T>& masses,
		       std::vector<LorentzVector<T>>& finalState, Bool_v multHi);

  T GenerateCosTheta(Int_v ptype, T pmod) const;

  // SPECIAL:  Generate N-body phase space using Kopylov algorithm
  void FillUsingKopylov(T initialMass,
			const std::vector<T>& masses,
			std::vector<LorentzVector<T>>& finalState);

  T BetaKopylov(size_t K) const;	// Copied from G4HadPhaseSpaceKopylov

private:
  const GXVMultiBodyMomDst* momDist;	// Buffers for selected distributions
  const G4VTwoBodyAngDst** angDist;	// Will be NULL for 3+body channels

  std::vector<Index_v<T>> kinds;	// Copy of particle_kinds list
  Index_v<T> vmultipl;		        // Final state size, for convenience
  int multiplicity;
  T bullet_ekin;			// Kinematics needed for distributions
  GXLorentzConvertor<T> toSCM;		// Handles complex rotations/transforms

  std::vector<T> modules;	// Buffers for generating momenta
  GXThreeVector<T> mom;

  static const T maxCosTheta;	// Cut for valid polar angle generation
  static const T oneOverE;	// Numeric value of 1/e for calculations
  static const T small;		// Cut for momentum/kinematics
  static constexpr size_t fTsize = vecCore::VectorSize<T>();
};

// Cut-offs and iteration limits for generation

template <typename T>
const T GXCascadeFinalStateAlgorithm<T>::maxCosTheta = 0.9999;

template <typename T>
const T GXCascadeFinalStateAlgorithm<T>::oneOverE = 0.3678794;   

template <typename T>
const T GXCascadeFinalStateAlgorithm<T>::small = 1.e-10;

// Constructor and destructor

template <typename T>
GXCascadeFinalStateAlgorithm<T>::GXCascadeFinalStateAlgorithm()
  : GXVCascadeAlgorithm<T>("GXCascadeFinalStateAlgorithm")
  , momDist(0)
  , angDist(0)
  , vmultipl(0)
  , multiplicity(0)
  , bullet_ekin(0.)
{
  angDist = new const G4VTwoBodyAngDst*[fTsize];
}

template <typename T>
GXCascadeFinalStateAlgorithm<T>::~GXCascadeFinalStateAlgorithm() { }

template <typename T>
void GXCascadeFinalStateAlgorithm<T>::SetVerboseLevel(int verbose) {
  GXVCascadeAlgorithm<T>::SetVerboseLevel(verbose);
  //GXMultiBodyMomentumDist::setVerboseLevel(verbose);
  GXTwoBodyAngularDist<T>::setVerboseLevel(verbose);
  toSCM.setVerbose(verbose);
}

// Algorithm declarations

template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
Configure(GXInuclElementaryParticle<T>* bullet,
	  GXInuclElementaryParticle<T>* target,
	  const std::vector<Int_v>& particle_kinds) {
  if (GetVerboseLevel()>1)
    std::cerr << " >>> " << GetName() << "::Configure()\n";

  // Identify initial and final state (if two-body) for algorithm selection
  multiplicity = particle_kinds.size();
  Int_v is = bullet->type() * target->type();
  Int_v fs(0);
  vecCore::MaskedAssign(fs, vmultipl == Index_v<T>(2), particle_kinds[0] * particle_kinds[1]);

  ChooseGenerators(is, fs);

  // Save kinematics for use with distributions
  SaveKinematics(bullet, target);

  // Save particle types for use with distributions
  kinds = particle_kinds;
}


template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
ChooseGenerators(Int_v is, Int_v fs) {
  if (GetVerboseLevel()>1) 
    std::cerr << " >>> " << GetName() << "::ChooseGenerators"
	      << " is " << is << " fs " << fs <<"\n";

  // Get generators for momentum and angle
  if (G4CascadeParameters::usePhaseSpace()) momDist = 0;
  else momDist = GXMultiBodyMomentumDist::GetDist(is, vmultipl);

  for(size_t i = 0; i < fTsize; ++i) {
    int isi = Get(is,i);
    int fsi = Get(fs,i);
    if (fsi > 0 && multiplicity == 2) {
      int kw = (fsi == isi) ? 1 : 2;
      angDist[i] = GXTwoBodyAngularDist<T>::GetDist(isi, fsi, kw);
    } else if (multiplicity == 3) {
      angDist[i] = GXTwoBodyAngularDist<T>::GetDist(isi);
    } else {
      angDist[i] = 0;
    }
  }

  if (GetVerboseLevel()>1) {
    std::cerr << " momDist: " << (momDist ? momDist->GetName().c_str() : "") << " - angDist: [";
    for(size_t i = 0; i < fTsize; ++i) {
      if (i>0) std::cerr <<"; ";
      std::cerr << (angDist[i] ? angDist[i]->GetName().c_str() : "none") <<"\n";
    }
    std::cerr << "]\n";
  }
}

template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
void GXCascadeFinalStateAlgorithm<T>::
FillMagnitudes(T initialMass, const std::vector<T>& masses) {
  assert(isHomogeneous(vmultipl));
  if (GetVerboseLevel()>1) 
    std::cerr << " >>> " << GetName() << "::FillMagnitudes\n";

  modules.clear();		// Initialization and sanity checks
  if (!momDist) return;

  modules.resize(multiplicity, 0.);	// Pre-allocate to avoid resizing

  T mass_last = masses.back();
  T pmod = 0.;

  if (GetVerboseLevel() > 3){
    std::cerr << " knd_last " << kinds.back() << " mass_last " << mass_last << "\n";
  }

  size_t itry = -1;
  while (++itry < itry_max) {		/* Loop checking 08.06.2015 MHK */
    if (GetVerboseLevel() > 3) {
      std::cerr << " itry in fillMagnitudes " << itry <<"\n";
    }

    T eleft = initialMass;

    size_t i;	// For access outside of loop
    size_t maxmult = multiplicity;
    for (i=0; i < maxmult-1; ++i) {
      pmod = momDist->GetMomentum(kinds[i], bullet_ekin);

      if (pmod < small) break;
      eleft -= std::sqrt(pmod*pmod + masses[i]*masses[i]);

      if (GetVerboseLevel() > 3) {
	std::cerr << " kp " << kinds[i] << " pmod " << pmod
		  << " mass2 " << masses[i]*masses[i] << " eleft " << eleft
		  << "\n x1 " << eleft - mass_last <<"\n";
      }

      if (eleft <= mass_last) break;

      modules[i] = pmod;
    }

    if (i < maxmult-1) continue;	// Failed to generate full kinematics

    T plast = eleft * eleft - masses.back()*masses.back();
    if (GetVerboseLevel() > 2) std::cerr << " plast ** 2 " << plast <<"\n";

    if (plast <= small) continue;	// Not enough momentum left over

    plast = std::sqrt(plast);		// Final momentum is what's left over
    modules.back() = plast;

    if (multiplicity > 3 || satisfyTriangle(modules)) break;	// Successful
  }	// while (itry < itry_max)

  if (itry >= itry_max) {		// Too many attempts
    if (GetVerboseLevel() > 2) {
      std::cerr << " Unable to generate momenta for multiplicity "
	<< multiplicity <<"\n";
    }
    modules.clear();		// Something went wrong, throw away partial
  }
}


// For three-body states, check kinematics of momentum moduli

template <typename T>
vecCore::Mask_v<T> GXCascadeFinalStateAlgorithm<T>::
satisfyTriangle(const std::vector<T>& pmod) const {
  if (GetVerboseLevel() > 3) 
    std::cerr << " >>> " << GetName() << "::satisfyTriangle\n";

  return ( (pmod.size() != 3) ||
	   !(pmod[0] < math::Abs(pmod[1] - pmod[2]) ||
	     pmod[0] > pmod[1] + pmod[2] ||
	     pmod[1] < math::Abs(pmod[0] - pmod[2]) ||
	     pmod[1] > pmod[0] + pmod[2] ||
	     pmod[2] < math::Abs(pmod[0] - pmod[1]) ||
	     pmod[2] > pmod[1] + pmod[0])
	   );
}


// Generate momentum directions into final state

template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
FillDirections(T initialMass, const std::vector<T>& masses,
	       std::vector<LorentzVector<T>>& finalState)
{
  if (GetVerboseLevel()>1) 
    std::cerr << " >>> " << GetName() << "::FillDirections.\n";

  finalState.clear();			// Initialization and sanity check
  size_t maxmult = multiplicity;
  finalState.resize(maxmult);

  // test
  if (modules.size() != maxmult) {
    std::cerr << " >>> " << GetName() << "::FillDirections() error: modules.size() != maxmult!\n";
    return; // ??? FIXME
  }

  // Different order of processing for three vs. N body
  Bool_v mult3 = (multiplicity == 3);
  Bool_v multHi = (multiplicity > 3);
  if (!vecCore::MaskEmpty(mult3))
    FillDirThreeBody(initialMass, masses, finalState, mult3);
  if (!vecCore::MaskEmpty(multHi))
    FillDirManyBody(initialMass, masses, finalState, multHi);
}


template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
FillDirThreeBody(T initialMass, const std::vector<T>& masses,
		 std::vector<LorentzVector<T>>& finalState, Bool_v mult3)
{
  if (GetVerboseLevel() > 1) {
    std::cerr << " >>> " << GetName() << "::FillDirThreeBody\n";
  }

  T costh = GenerateCosTheta(kinds[2], modules[2]);
  // FIXME: should use a local LorVector here instead of finalstate[2]?
  vecCore::MaskedAssign(finalState[2], mult3, generateWithFixedTheta(costh, modules[2], masses[2]));
  vecCore::MaskedAssign(finalState[2], mult3, toSCM.rotate(finalState[2]));	// Align target axis

  // Generate direction of first particle
  costh = -0.5 * (modules[2]*modules[2] + modules[0]*modules[0] - modules[1]*modules[1]) / modules[2] / modules[0];

  if (std::fabs(costh) >= maxCosTheta) {  // Bad kinematics; abort generation
    finalState.clear();
    return;
  }

  // Report success
  if (GetVerboseLevel()>2) {
    std::cerr << " ok for mult 3\n";
  }

  // First particle is at fixed angle to recoil system
  finalState[0] = generateWithFixedTheta(costh, modules[0], masses[0]);
  finalState[0] = toSCM.rotate(finalState[2], finalState[0]);

  // Remaining particle is constrained to recoil from entire rest of system
  finalState[1].Set(0.,0.,0.,initialMass);
  finalState[1] -= finalState[0] + finalState[2];
}

template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
FillDirManyBody(T initialMass, const std::vector<T>& masses,
		std::vector<LorentzVector<T>>& finalState, Bool_v multHi)
{
  if (GetVerboseLevel()>1) {
    std::cerr << " >>> " << GetName() << "::FillDirManyBody\n";
  }

  // Fill all but the last two particles randomly
  T costh = 0.;

  finalState.resize(multiplicity);

  for (size_t i = 0; i < multiplicity-2; ++i) {
    costh = GenerateCosTheta(kinds[i], modules[i]);
    finalState[i] = generateWithFixedTheta(costh, modules[i], masses[i]);
    finalState[i] = toSCM.rotate(finalState[i]);	// Align target axis
  }

  // Total momentum so far, to compute recoil of last two particles
  LorentzVector<T> psum =
    std::accumulate(finalState.begin(), finalState.end()-2, LorentzVector<T>());
  T pmod = psum.Vect().Mag();

  costh = -0.5 * (pmod*pmod +
		  modules[multiplicity-2]*modules[multiplicity-2] -
		  modules[multiplicity-1]*modules[multiplicity-1])
    / pmod / modules[multiplicity-2];

  if (GetVerboseLevel() > 2) {
    std::cerr <<  " ct last " << costh <<"\n";
  }

  if (std::fabs(costh) >= maxCosTheta) {  // Bad kinematics; abort generation
    finalState.clear();
    return;
  }

  // Report success
  if (GetVerboseLevel()>2) std::cerr << " ok for mult " << multiplicity <<"\n";

  // First particle is at fixed angle to recoil system
  finalState[multiplicity-2] =
    generateWithFixedTheta(costh, modules[multiplicity-2],
			   masses[multiplicity-2]);
  finalState[multiplicity-2] = toSCM.rotate(psum, finalState[multiplicity-2]);

  // Remaining particle is constrained to recoil from entire rest of system
  finalState[multiplicity-1].Set(0.,0.,0.,initialMass);
  finalState[multiplicity-1] -= psum + finalState[multiplicity-2];
}


// Generate polar angle for three- and multi-body systems

template <typename T>
T GXCascadeFinalStateAlgorithm<T>::
GenerateCosTheta(Int_v ptype, T pmod) const {
  constexpr size_t vsize = vecCore::VectorSize<T>();
  using Int_T = typename vecCore::backend::VcSimdArray<vsize>::Int_v;
  GXPowVec<T,Int_T>* ppow = GXPowVec<T, Int_T>::GetInstance();

  if (GetVerboseLevel() > 2) {
    std::cerr << " >>> " << GetName() << "::GenerateCosTheta " << ptype
	      << " " << pmod <<"\n";
  }

  if (multiplicity == 3) {		// Use distribution for three-body
    T result(0.);
    for(size_t i = 0; i < fTsize; ++i) {
      double ekin = Get(bullet_ekin, i);
      int itype = Get(ptype, i);
      double costhi = angDist[i]->GetCosTheta(ekin, itype);
      Set(result, i, costhi);
    }
    return result;
  }

  //.. Throw multi-body distribution
  // original: G4double p0 = ptype<3 ? 0.36 : 0.25;	// Nucleon vs. everything else
  T p0(0.36);
  for(size_t i = 0; i < vsize; ++i) {
    if(Get(ptype,i) >= 3) Set(p0, i, 0.25);
  }

  T alf = T(1.0) / (p0 * (p0 - (pmod + p0) * ppow->ExpAVec(-pmod / p0)));
  T sinth(2.0);

  int itry1 = 0;		/* Loop checking 08.06.2015 MHK */
  Bool_v done(false);
  do {
    T s1 = pmod * inuclRndm<T>();
    T s2 = alf * oneOverE * p0 * inuclRndm<T>();
    T salf = s1 * alf * ppow->ExpAVec(-s1 / p0);
    if (GetVerboseLevel() > 3) {
      std::cerr << " s1 * alf * GXExp(-s1 / p0) " << salf << " s2 " << s2 <<"\n";
    }

    //if (salf > s2) sinth = s1 / pmod;
    vecCore::MaskedAssign(sinth, !done && salf > s2, s1 / pmod);
    done = math::Abs(sinth) <= maxCosTheta;
  } while (!vecCore::MaskFull(done) && ++itry1 < itry_max);

  //if (GetVerboseLevel() > 3)
  std::cerr << "GXCascadeFSAlgorithm: GenerateCosTheta(): itry1=" << itry1 << " sinth=" << sinth <<"\n";

  if (itry1 == itry_max) {
    if (GetVerboseLevel() > 2) {
      std::cerr << " high energy angles generation: itry1 " << itry1 <<"\n";;
    }
    vecCore::MaskedAssign(sinth, !done, 0.5 * inuclRndm<T>());
  }

  // Convert generated sin(theta) to cos(theta) with random sign
  T costh = math::Sqrt(1.0 - sinth * sinth);
  vecCore::MaskedAssign(costh, inuclRndm<T>() > 0.5, -costh);

  return costh;
}


template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
GenerateTwoBody(T initialMass, std::vector<T> const& masses,
		std::vector<LorentzVector<T>>& finalState)
{
  if (GetVerboseLevel() > 1) 
    std::cerr << " >>> " << GetName() << "::GenerateTwoBody.\n";

  finalState.clear();		// Initialization and sanity checks

  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskFull(multiplicity != Index_v<T>(2))) return;
  }

  {
    // Generate momentum vector in CMS for back-to-back particles
    const T pscm = this->TwoBodyMomentum(initialMass, masses[0], masses[1]);

    T costh(2.);
    for(size_t i = 0; i < fTsize; ++i) {
      // GL.. to be fully vectorized eventually
      double costhi = angDist[i] ? angDist[i]->GetCosTheta(Get(bullet_ekin,i), Get(pscm,i)) : (2. * inuclRndm<double>() - 1.);
      Set(costh, i, costhi);
    }

    //mom.setRThetaPhi(pscm, Acos(costh), this->UniformPhi());
    mom.SetMagCosThPhi(pscm, costh, this->UniformPhi());

    if (GetVerboseLevel() > 3) {		// Copied from old G4EPCollider
      std::cerr << " Particle kinds = " << kinds[0] << " , " << kinds[1]
		<< "\n pmod " << pscm
		<< "\n before rotation px " << mom.x() << " py " << mom.y()
		<< " pz " << mom.z() <<"\n";
    }
  }

  finalState.resize(2);				// Allows filling by index

  finalState[0].SetVectMag(mom, masses[0]);
  finalState[0] = toSCM.rotate(finalState[0]);

  if (GetVerboseLevel()>3) {		// Copied from old G4EPCollider
    std::cerr << " after rotation px " << finalState[0].x() << " py "
	      << finalState[0].y() << " pz " << finalState[0].z() <<"\n";;
  }

  finalState[1].SetVectMag(-finalState[0].Vect(), masses[1]);
}

// N-body generation uses momentum-modulus distribution, computed angles

template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
GenerateMultiBody(T initialMass, const std::vector<T>& masses,
		  std::vector<LorentzVector<T>>& finalState) {

  if (GetVerboseLevel()>1)
    std::cerr << " >>> " << GetName() << "::GenerateMultiBody\n";

  if (G4CascadeParameters::usePhaseSpace()) {
    FillUsingKopylov(initialMass, masses, finalState);
    return;
  }

  finalState.clear();		// Initialization and sanity checks

  Bool_v done = (multiplicity < Index_v<T>(3));
  if (vecCore::EarlyReturnAllowed()) {
    if (vecCore::MaskFull(done) || momDist == NULL) return;
  }

  int itry = -1;		/* Loop checking 08.06.2015 MHK */
  while (!vecCore::MaskFull(done) && ++itry < itry_max) {
    FillMagnitudes(initialMass, masses);
    FillDirections(initialMass, masses, finalState);
    done |= (T(finalState.size()) == multiplicity);
  }
}


template <typename T>
void GXCascadeFinalStateAlgorithm<T>::
FillUsingKopylov(T initialMass,
		 std::vector<T> const& masses,
		 std::vector<LorentzVector<T>>& finalState)
{
  if (GetVerboseLevel()>2)
    std::cerr << " >>> " << GetName() << "::FillUsingKopylov.\n";

  finalState.clear();

  size_t N = masses.size();
  finalState.resize(N);

  T zero(0.0);
  T mtot = std::accumulate(masses.begin(), masses.end(), zero);
  T mu = mtot;
  T Mass = initialMass;
  T temp = Mass - mtot;
  T recoilMass;
  GXThreeVector<T> momV, boostV;		// Buffers to reduce memory churn
  LorentzVector<T> recoil(0.0, 0.0, 0.0, Mass);

  for (size_t k = N-1; k>0; --k) {
    mu -= masses[k];
    temp *= (k > 1) ? BetaKopylov(k) : zero;

    recoilMass = mu + temp;

    boostV = recoil.BoostVector();	// Previous system's rest frame

    {
      // Create momentum with a random direction isotropically distributed
      // FIXME (TKoi|MKelsey):  Should theta distribution use Bertini fit function?

      //momV.setRThetaPhi(TwoBodyMomentum(Mass,masses[k],recoilMass), UniformTheta(), UniformPhi());
      const T totmom = this->TwoBodyMomentum(Mass, masses[k], recoilMass);
      const T costh = (2. * inuclRndm<T>() - 1.);
      momV.SetMagCosThPhi( totmom, costh, this->UniformPhi());
    }

    finalState[k].SetVectMag(momV, masses[k]);
    recoil.SetVectMag(-momV, recoilMass);

    finalState[k].Boost(boostV);
    recoil.Boost(boostV);
    Mass = recoilMass;
  }
  
  finalState[0] = recoil;
}

template <typename T>
T GXCascadeFinalStateAlgorithm<T>::
BetaKopylov(size_t K) const
{
  constexpr size_t vsize = vecCore::VectorSize<T>();
  using Int_v = typename vecCore::backend::VcSimdArray<vsize>::Int_v;
  GXPowVec<T,Int_v>* g4pow = GXPowVec<T,Int_v>::GetInstance();

  Index_v<T> N(3 * K - 5);
  const T xN( N );
  const T one(1.0);

  // original: Fmax = std::sqrt(g4pow->powN(xN/(xN+1.),N)/(xN+1.)); 
  const T oneOverXnPlus1 = one / (xN + one);
  const T Fmax = math::Sqrt(g4pow->PowN(xN * oneOverXnPlus1, N) * oneOverXnPlus1); 

  // Loop checking 08.06.2015 MHK
  T F, chi;
  vecCore::Mask_v<T> done(false);
  do {
    vecCore::MaskedAssign( chi, !done, inuclRndm<T>());
    vecCore::MaskedAssign( F, !done, math::Sqrt(g4pow->PowN(chi,N) * (one - chi)));
    done = (Fmax * inuclRndm<T>() < F);
  } while ( !vecCore::MaskFull(done) );
  return chi;
}

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXCascadeFinalStateAlgorithm_hh
