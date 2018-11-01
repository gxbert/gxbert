//
// @File: GXCascadeGenerator.hh
//
#ifndef GXCascadeGenerator_HH
#define GXCascadeGenerator_HH 1

#include "LorentzVector.hh"
#include "GXHadPhaseSpaceKopylov.hh"
#include "GXHadPhaseSpaceGenbod.hh"
#include "GXHadPhaseSpaceNBodyAsai.hh"
#include <vector>

class GXParticleDefinition;

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

template <typename T> class GXVCascadeAlgorithm;

template <typename T>
class GXCascadeGenerator {

public:
  using Bool_v = vecCore::Mask_v<T>;

  // Flags to select algorithm by code in constructor
  enum Algorithm { NONE=0, Kopylov=1, GENBOD=2, NBody=3 };

  // Specify "standard" algorithm by code or by object (takes ownership)
  GXCascadeGenerator(Algorithm alg = Kopylov, int verbose = 0);
  GXCascadeGenerator(GXVCascadeAlgorithm<T>* alg, int verbose = 0);
  virtual ~GXCascadeGenerator();

  // Enable (or disable if 0) diagnostic messages
  void SetVerboseLevel(int verbose);

  const std::string& GetAlgorithmName() const;

  // Initial state (rest mass) and list of final masses
  bool Generate(T initialMass,
                std::vector<T> const& masses,
                std::vector<LorentzVector<T>>& finalState);

  // Initial state particle and list of final masses
  bool Generate(GXParticleDefinition const* initialPD,
                std::vector<T> const& masses,
      	        std::vector<LorentzVector<T>>& finalState);

  // Initial state (frame) and list of final masses
  // Final state particles will be boosted to initial-state frame
  bool Generate(LorentzVector<T> const& initialState,
                std::vector<T> const& masses,
                std::vector<LorentzVector<T>>& finalState);

protected:
  // Special case for one-body final state
  bool GenerateOneBody(T initialMass,
                       const std::vector<T>& masses,
                       std::vector<LorentzVector<T>>& finalState) const;

  // Special function used by constructor for unrecognized algorithm code
  void ReportInvalidAlgorithm(Algorithm alg) const;
  void ReportMissingAlgorithm() const;

protected:
  // SPECIAL FUNCTION FOR SUBCLASSES: A subclass may implement a
  // collection of algorithms, to be switched on an event-by-event
  // basis.  This function allows the subclass to switch the "active"
  // algorithm before Generate() is called.
  //
  // If this function is used by the subclass, then the subclass has
  // ownership of _all_ instantiated algorithms, and should delete
  // them in its own dtor.  The subclass dtor must also call
  // UseAlgorithm(0) to set the base algorithm to a null pointer, to
  // prevent a double-delete error.
  void UseAlgorithm(GXVCascadeAlgorithm<T>* alg) { theAlgorithm = alg; }

  int verboseLevel;
  GXVCascadeAlgorithm<T>* theAlgorithm;
};

template <typename T>
GXCascadeGenerator<T>::GXCascadeGenerator(Algorithm alg, G4int verbose)
  : verboseLevel(verbose), theAlgorithm(0) {
  switch (alg) {
  case Kopylov: theAlgorithm = new GXHadPhaseSpaceKopylov(verboseLevel); break;
  case GENBOD: theAlgorithm = new GXHadPhaseSpaceGenbod(verboseLevel); break;
  case NBody: theAlgorithm = new GXHadPhaseSpaceNBodyAsai(verboseLevel); break;
  case NONE: theAlgorithm = 0; break;	// User may explicitly set no algorithm
  default: ReportInvalidAlgorithm(alg);
  }

  if (verboseLevel) {
    G4cout << " >>> GXHadDecayGenerator";
    if (theAlgorithm) G4cout << " using " << theAlgorithm->GetName();
    G4cout << G4endl;
  }
}

template <typename T>
GXCascadeGenerator<T>::GXCascadeGenerator(GXVCascadeAlgorithm<T>* alg,
					  int verbose)
  : verboseLevel(verbose), theAlgorithm(alg) {
  if (verboseLevel>1) {
    std::cerr << " >>> GXCascadeGenerator";
    if (theAlgorithm) std::cerr << " using " << theAlgorithm->GetName();
    std::cerr <<"\n";
  }
}

template <typename T>
GXCascadeGenerator<T>::~GXCascadeGenerator()
{
  if (theAlgorithm) {
    delete theAlgorithm;
    theAlgorithm = 0;
  }
}

template <typename T>
void GXCascadeGenerator<T>::SetVerboseLevel(int verbose)
{
  verboseLevel = verbose;
  if (verboseLevel>1)
    std::cerr << " >>> GXCascadeGenerator::SetVerboseLevel("<< verbose <<")\n";
  if (theAlgorithm) theAlgorithm->SetVerboseLevel(verbose);
}

// Initial state (rest mass) and list of final masses

template <typename T>
bool GXCascadeGenerator<T>::Generate(T initialMass,
				     const std::vector<T>& masses,
				     std::vector<LorentzVector<T>>& finalState)
{
  if (verboseLevel>1)
    std::cerr << " >>> GXCascadeGenerator::Generate (mass) [EMPTY!!]\n";

  if (!theAlgorithm) ReportMissingAlgorithm();

  if (masses.size() == 1U)
    return GenerateOneBody(initialMass, masses, finalState);

  theAlgorithm->Generate(initialMass, masses, finalState);

  return !finalState.empty();		// Generator failure returns empty state
}

// Initial state particle and list of final masses

template <typename T>
bool GXCascadeGenerator<T>::
Generate(GXParticleDefinition const* initialPD,
	 std::vector<T> const& masses,
	 std::vector<LorentzVector<T>>& finalState)
{
  if (verboseLevel>1) std::cerr << " >>> GXCascadeGenerator::Generate (particle).\n";
  return (initialPD && Generate(T(initialPD->GetPDGMass()), masses, finalState));
}

// Final state particles will be boosted to initial-state frame

template <typename T>
bool GXCascadeGenerator<T>::
Generate(LorentzVector<T> const& initialState,
	 std::vector<T> const& masses,
	 std::vector<LorentzVector<T>>& finalState)
{
  if (verboseLevel>1) {
    std::cerr << " >>> GXCascadeGenerator::Generate (frame).\n";
  }
  bool good = Generate(initialState.m(), masses, finalState);
  if (good) {
    GXThreeVector<T> bv = initialState.boostVector();
    for (size_t i = 0; i < finalState.size(); ++i) {
      finalState[i].boost(bv);
    }
  }

  return good;
}

// Handle special case of "one body decay" (used for kaon mixing)

template <typename T>
bool GXCascadeGenerator<T>::
GenerateOneBody(T initialMass,
		std::vector<T> const& masses,
		std::vector<LorentzVector<T>>& finalState) const
{
  if (verboseLevel>1) std::cerr << " >>> GXCascadeGenerator<T>::GenerateOneBody...\n";

  // Initialization and sanity checks
  finalState.clear();

  if (masses.size() != 1U) return false;	// Should not have been called
  vecCore::Mask_v<T> done = ( math::Abs(initialMass - masses[0]) > T(CLHEP::eV) );
  if(vecCore::MaskFull(done)) return false;
  assert(!vecCore::MaskEmpty(done) &&
	 ">>> GXCascadeGenerator<T>::GenerateOneBody(): Warning: this line should never be reached!\n");

  if (verboseLevel>2) std::cerr << " finalState mass = " << masses[0] <<"\n";

  finalState.push_back(LorentzVector<T>(0.,0.,0.,masses[0]));
  return true;
}

template <typename T>
void GXCascadeGenerator<T>::
ReportMissingAlgorithm() const
{
  if (verboseLevel>1) {
    std::cerr << "GXCascadeGenerator<T>: no algorithm specified.\n";
  }
  throw GXHadronicException(__FILE__, __LINE__, "Null algorithm pointer");
}

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXCascadeGenerator_HH
