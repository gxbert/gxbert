//
// @File: GXCascadeGenerator.hh
//
#ifndef GXCascadeGenerator_HH
#define GXCascadeGenerator_HH 1

//#include "globals.hh"
#include "LorentzVector.hh"
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
  Bool_v Generate(T initialMass,
		  const std::vector<T>& masses,
		  std::vector<LorentzVector<T>>& finalState);

  // Initial state particle and list of final masses
  Bool_v Generate(const GXParticleDefinition* initialPD,
		  const std::vector<T>& masses,
		  std::vector<LorentzVector<T>>& finalState);

  // Initial state (frame) and list of final masses
  // Final state particles will be boosted to initial-state frame
  Bool_v Generate(const LorentzVector<T>& initialState,
		  const std::vector<T>& masses,
		  std::vector<LorentzVector<T>>& finalState);

protected:
  // Special case for one-body final state
  Bool_v GenerateOneBody(T initialMass,
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
GXCascadeGenerator<T>::GXCascadeGenerator(GXVCascadeAlgorithm<T>* alg,
					  int verbose)
  : verboseLevel(verbose), theAlgorithm(alg) {
  if (verboseLevel) {
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
  if (theAlgorithm) theAlgorithm->SetVerboseLevel(verbose);
}

template <typename T>
vecCore::Mask_v<T> GXCascadeGenerator<T>::Generate(T initialMass,
						   const std::vector<T>& masses,
						   std::vector<LorentzVector<T>>& finalState)
{
  if (verboseLevel) 
    G4cout << " >>> GXCascadeGenerator::Generate (mass) [EMPTY!!]" << G4endl;
  /*
  if (!theAlgorithm) ReportMissingAlgorithm();

  if (masses.size() == 1U)
    return GenerateOneBody(initialMass, masses, finalState);

  theAlgorithm->Generate(initialMass, masses, finalState);
  */
  return !finalState.empty();		// Generator failure returns empty state
}

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif	/* GXCascadeGenerator_HH */
