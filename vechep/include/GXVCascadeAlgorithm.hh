//
// @File: GXVCascadeAlgorithm.hh
//
// 20180713  Guilherme Lima -- Created, based on T.Koi's GXVHadDecayAlgorithm

#ifndef GXVCascadeAlgorithm_HH
#define GXVCascadeAlgorithm_HH 1

#include "GXInuclSpecialFunctions.hh"
#include "LorentzVector.hh"
#include "GXThreeVector.hh"

#include <vector>
#include <string>
#include <numeric>

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

using namespace GXInuclSpecialFunctions;

template <typename T>
class GXVCascadeAlgorithm {

  using Bool_v = vecCore::Mask_v<T>;

public:
  GXVCascadeAlgorithm(const std::string& algName, int verbose=0)
    : name(algName), verboseLevel(verbose) {;}
  virtual ~GXVCascadeAlgorithm() {;}

  // Initial state (rest mass) and list of final masses
  void Generate(T initialMass,
		std::vector<T> const& masses,
		std::vector<LorentzVector<T>>& finalState);

  // Enable (or disable if 0) diagnostic messages (subclass may overload)
  void  SetVerboseLevel(int verbose)
  {
    verboseLevel = verbose;
    if (GetVerboseLevel() > 2)
      std::cerr << " >>> "<< GetName() << "::SetVerboseLevel("<< verbose <<")\n";
  }

  int GetVerboseLevel() const { return verboseLevel; }
  const std::string& GetName() const { return name; }
  
protected:
  // Subclasses MUST implement these functions
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  virtual void GenerateTwoBody(T initialMass,
			       const std::vector<T>& masses,
			       std::vector<LorentzVector<T>>& finalState) = 0;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  virtual void GenerateMultiBody(T initialMass,
				 const std::vector<T>& masses,
				 std::vector<LorentzVector<T>>& finalState) = 0;

  // Validate kinematics (e.g., limit number of final state particles)
  // Subclasses may override or call back to this function
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  virtual Bool_v IsDecayAllowed(T initialMass,
				std::vector<T> const& masses) const;

  // Two-body momentum function (c.f. PDK from CERNLIB W505)
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T TwoBodyMomentum(T M0, T M1, T M2) const;

  // Convenience functions for uniform angular distributions
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T UniformCosTheta() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T UniformTheta() const;

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T UniformPhi() const;

  // Utility to dump vector contents to line of output
  void PrintVector(const std::vector<T>& v,
		   const std::string& name,
		   std::ostream& os) const;

private:
  std::string name;
  int verboseLevel;
};

  // Base class does very simple validation of configuration

  template <typename T>
  vecCore::Mask_v<T> GXVCascadeAlgorithm<T>::
  IsDecayAllowed(T initialMass, std::vector<T> const& masses) const
  {
    Bool_v okay( masses.size() >= 2 );
    okay = okay & (initialMass > T(0.)) & (initialMass >= std::accumulate(masses.begin(), masses.end(), T(0.)));

    if (verboseLevel>1) {
      std::cerr << GetName() << "::IsDecayAllowed? initialMass " << initialMass
		<< " - " << masses.size() << " masses sum "
		<< std::accumulate(masses.begin(), masses.end(), T(0.)) <<"\n";

      if (verboseLevel>1) PrintVector(masses, " ", std::cerr);

      std::cerr << " Returning " << okay <<"\n";
    }

    return okay;
  }


  // Momentum function (c.f. PDK() function from CERNLIB W515)
  template <typename T>
  T GXVCascadeAlgorithm<T>::
  TwoBodyMomentum(T M0, T M1, T M2) const
  {
    using CLHEP::GeV;
    using CLHEP::MeV;

    T PSQ = (M0 + M1 + M2) * (M0 + M1 - M2) * (M0 - M1 + M2) * (M0 - M1 - M2);
    if (!vecCore::MaskEmpty(PSQ < 0.)) {
      std::cerr << GetName() << ":  problem of decay of M(GeV) " << M0/GeV
		<< " to M1(GeV) " << M1/GeV << " and M2(GeV) " << M2/GeV
		<< " PSQ(MeV) " << PSQ/MeV << " < 0\n";

      // exception only if the problem is numerically significant
      if (!vecCore::MaskEmpty(PSQ < -CLHEP::eV)) {
	throw GXHadronicException(__FILE__, __LINE__,"Error in decay kinematics");
      }

      vecCore::MaskedAssign(PSQ, PSQ < 0., T(0.));
    }

    T result = math::Sqrt(PSQ) / (2. * M0);
    if (verboseLevel > 1) std::cerr<<"TwoBodyMom(): returning T="<< result <<"\n";
    return result;
  }

  // Convenience functions for uniform angular distributions

  template <typename T>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T GXVCascadeAlgorithm<T>::UniformCosTheta() const
  {
    return 2.0 * inuclRndm<T>() - 1.0;
  }

  // Note: keeping original name for backward compatibility, but distribution is not flat in theta! (and Acos() is slow)
  template <typename T>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T GXVCascadeAlgorithm<T>::UniformTheta() const
  {
    //return 0.5 * randomPHI();  // this is flat in theta!
    return std::acos(2.0 * inuclRndm<T>() - 1.0);
  }

  template <typename T>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T GXVCascadeAlgorithm<T>::UniformPhi() const
  {
    return randomPHI<T>();
    //return CLHEP::twopi * G4UniformRand();
  }


  // Dump contents of vector to output

  template <typename T>
  void GXVCascadeAlgorithm<T>::
  PrintVector(const std::vector<T>& v,
	      const std::string& vname, std::ostream& os) const
  {
    os << " " << vname << "(" << v.size() << ") ";
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    os <<"\n";
  }

  // Initial state (rest mass) and list of final masses
  template <typename T>
  void GXVCascadeAlgorithm<T>::
  Generate(T initialMass,
	   std::vector<T> const& masses,
	   std::vector<LorentzVector<T>>& finalState)
  {
    if (verboseLevel>1) std::cerr << GetName() << "::Generate.\n";

    // Initialization and sanity check
    finalState.clear();
    Bool_v doit = IsDecayAllowed(initialMass, masses);
    if(vecCore::MaskEmpty(doit)) return;
    assert(vecCore::MaskFull(doit) && "GXVCascadeAlgorithm(): mixed decay reasonability in lanes!");

    // Allow different procedures for two-body or N-body distributions
    if (masses.size() == 2U) {
      GenerateTwoBody(initialMass, masses, finalState);
    }
    else {
      GenerateMultiBody(initialMass, masses, finalState);
    }
  }

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif	/* GXVCascadeAlgorithm_HH */
