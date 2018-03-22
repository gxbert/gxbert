//
// Geant4 header GXFragment
//
// by V. Lara (May 1998)
//
// Modifications:
// 03.05.2010 V.Ivanchenko General cleanup of inline functions: objects 
//            are accessed by reference; remove double return 
//            tolerance of excitation energy at modent it is computed;
//            safe computation of excitation for exotic fragments
// 18.05.2010 V.Ivanchenko added member theGroundStateMass and inline
//            method which allowing to compute this value once and use 
//            many times
// 26.09.2010 V.Ivanchenko added number of protons, neutrons, proton holes
//            and neutron holes as members of the class and Get/Set methods;
//            removed not needed 'const'; removed old debug staff and unused
//            private methods; add comments and reorder methods for 
//            better reading

#ifndef GXFragment_h
#define GXFragment_h 1

#include "globals.hh"
#include "GXAllocator.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "GXNucleiProperties.hh"
//#include "G4PhysicalConstants.hh"
//#include "CLHEP/Units/PhysicalConstants.h"
//#include "G4Proton.hh"
//#include "G4Neutron.hh"
#include <vector>

//class G4ParticleDefinition;

class GXFragment;     
typedef std::vector<GXFragment*> GXFragmentVector;

class GXFragment 
{
public:

  // ============= CONSTRUCTORS ==================

  // Default constructor - obsolete
  GXFragment();

  // Destructor
  ~GXFragment();

  // Copy constructor
  GXFragment(const GXFragment &right);

  // A,Z and 4-momentum - main constructor for fragment
  GXFragment(G4int A, G4int Z, const G4LorentzVector& aMomentum);

  // 4-momentum and pointer to G4particleDefinition (for gammas, e-)
//  GXFragment(const G4LorentzVector& aMomentum, 
//	     const G4ParticleDefinition* aParticleDefinition);

  // ============= OPERATORS ==================
    
  GXFragment & operator=(const GXFragment &right);
  G4bool operator==(const GXFragment &right) const;
  G4bool operator!=(const GXFragment &right) const;

  friend std::ostream& operator<<(std::ostream&, const GXFragment&);

  //  new/delete operators are overloded to use GXAllocator
  inline void *operator new(size_t);
  inline void operator delete(void *aFragment);

  // ============= GENERAL METHODS ==================

  inline G4int GetZ_asInt() const;
  inline G4int GetA_asInt() const;
  inline void SetZandA_asInt(G4int Znew, G4int Anew);
  
  inline G4double GetExcitationEnergy() const;
  inline void SetExcEnergyAndMomentum(G4double eexc, const G4LorentzVector&);

  inline G4double GetGroundStateMass() const;
   
//  inline G4double GetBindingEnergy() const;
  
  inline const G4LorentzVector& GetMomentum() const;
  inline void SetMomentum(const G4LorentzVector& value);
  
  // computation of mass for any Z and A
  inline G4double ComputeGroundStateMass(G4int Z, G4int A) const;

  // extra methods
//  inline G4double GetSpin() const;
//  inline void SetSpin(G4double value);

//  inline G4int GetCreatorModelType() const;
//  inline void SetCreatorModelType(G4int value);

  // obsolete methods
  
//  inline G4double GetZ() const;
//  inline G4double GetA() const;
//  inline void SetZ(G4double value);
//  inline void SetA(G4double value);
  
  // ============= METHODS FOR PRE-COMPOUND MODEL ===============

  inline G4int GetNumberOfExcitons() const;
  
  inline G4int GetNumberOfParticles() const;
  inline G4int GetNumberOfCharged() const;
  inline void SetNumberOfExcitedParticle(G4int valueTot, G4int valueP);

  inline G4int GetNumberOfHoles() const;
  inline G4int GetNumberOfChargedHoles() const;
  inline void SetNumberOfHoles(G4int valueTot, G4int valueP=0);
  
  // these methods will be removed in future
  inline void SetNumberOfParticles(G4int value);
  inline void SetNumberOfCharged(G4int value);

  // ============= METHODS FOR PHOTON EVAPORATION ===============

//  inline G4int GetNumberOfElectrons() const;
//  inline void SetNumberOfElectrons(G4int value);

//  inline G4int GetFloatingLevelNumber() const;
//  inline void SetFloatingLevelNumber(G4int value);

//  inline const G4ParticleDefinition * GetParticleDefinition() const;
//  inline void SetParticleDefinition(const G4ParticleDefinition * p);

//  inline G4double GetCreationTime() const;
//  inline void SetCreationTime(G4double time);

  // GXFragment class is not responsible for creation and delition of 
  // G4NuclearPolarization object
//  inline G4NuclearPolarization* NuclearPolarization();
//  inline G4NuclearPolarization* GetNuclearPolarization() const;
//  inline void SetNuclearPolarization(G4NuclearPolarization*);

//  void SetAngularMomentum(const G4ThreeVector&);
//  G4ThreeVector GetAngularMomentum() const;

  // ============= PRIVATE METHODS ==============================

private:

  void ExcitationEnergyWarning();

  void NumberOfExitationWarning(const G4String&);

  inline void CalculateExcitationEnergy();

  inline void CalculateGroundStateMass();

  // ============= DATA MEMBERS ==================

  G4int theA;
  
  G4int theZ;
  
  G4double theExcitationEnergy;

  G4double theGroundStateMass;

  G4LorentzVector theMomentum;
  
  // Nuclear polarisation by default is nullptr
//  G4NuclearPolarization* thePolarization;

  // creator model type
//  G4int creatorModel;

  // Exciton model data members  
  G4int numberOfParticles;  
  G4int numberOfCharged;
  G4int numberOfHoles;
  G4int numberOfChargedHoles;

  // Gamma evaporation data members
//  G4int numberOfShellElectrons;
//  G4int xLevel;

//  const G4ParticleDefinition* theParticleDefinition;
  
//  G4double spin;
//  G4double theCreationTime;

  static const G4double minFragExcitation;
};

// ============= INLINE METHOD IMPLEMENTATIONS ===================

#if defined G4HADRONIC_ALLOC_EXPORT
  extern G4DLLEXPORT G4ThreadLocal GXAllocator<GXFragment> *pGXFragmentAllocator;
#else
  extern G4DLLIMPORT G4ThreadLocal GXAllocator<GXFragment> *pGXFragmentAllocator;
#endif

inline void * GXFragment::operator new(size_t)
{
  if (!pGXFragmentAllocator) { pGXFragmentAllocator = new GXAllocator<GXFragment>; }
  return (void*) pGXFragmentAllocator->MallocSingle();
}

inline void GXFragment::operator delete(void * aFragment)
{
  pGXFragmentAllocator->FreeSingle((GXFragment *) aFragment);
}

inline void GXFragment::CalculateExcitationEnergy()
{
  theExcitationEnergy = theMomentum.mag() - theGroundStateMass;
  if(theExcitationEnergy < minFragExcitation) { 
    if(theExcitationEnergy < -minFragExcitation) {  ExcitationEnergyWarning(); }
    theExcitationEnergy = 0.0;
  }
}

inline G4double 
GXFragment::ComputeGroundStateMass(G4int Z, G4int A) const
{
  return GXNucleiProperties::GetNuclearMass(A, Z); 
}
	 
inline void GXFragment::CalculateGroundStateMass() 
{
  theGroundStateMass = GXNucleiProperties::GetNuclearMass(theA, theZ);
}

inline G4int GXFragment::GetA_asInt() const
{
  return theA;
}

inline G4int GXFragment::GetZ_asInt()  const
{
  return theZ;
}

inline void GXFragment::SetZandA_asInt(G4int Znew, G4int Anew)
{
  theZ = Znew;
  theA = Anew;
  CalculateGroundStateMass();
}

inline G4double GXFragment::GetExcitationEnergy()  const
{
  return theExcitationEnergy;
}

inline G4double GXFragment::GetGroundStateMass() const
{
  return theGroundStateMass; 
}

inline void GXFragment::SetExcEnergyAndMomentum(G4double eexc, 
						const G4LorentzVector& v)
{
  theExcitationEnergy = eexc;
  theMomentum.set(0.0, 0.0, 0.0, theGroundStateMass + eexc);
  theMomentum.boost(v.boostVector());
}

/*
inline G4double GXFragment::GetBindingEnergy() const
{
  return (theA-theZ)*CLHEP::neutron_mass_c2 + theZ*CLHEP::proton_mass_c2 
    - theGroundStateMass;
}
*/

inline const G4LorentzVector& GXFragment::GetMomentum()  const
{
  return theMomentum;
}

inline void GXFragment::SetMomentum(const G4LorentzVector& value)
{
  theMomentum = value;
  CalculateExcitationEnergy();
}

/*
inline G4double GXFragment::GetZ()  const
{
  return G4double(theZ);
}

inline G4double GXFragment::GetA() const
{
  return G4double(theA);
}

inline void GXFragment::SetZ(const G4double value)
{
  theZ = G4lrint(value);
  CalculateGroundStateMass();
}

inline void GXFragment::SetA(const G4double value)
{
  theA = G4lrint(value);
  CalculateGroundStateMass();
}
*/

inline G4int GXFragment::GetNumberOfExcitons()  const
{
  return numberOfParticles + numberOfHoles;
}

inline G4int GXFragment::GetNumberOfParticles()  const
{
  return numberOfParticles;
}

inline G4int GXFragment::GetNumberOfCharged()  const
{
  return numberOfCharged;
}

inline 
void GXFragment::SetNumberOfExcitedParticle(G4int valueTot, G4int valueP)
{
  numberOfParticles = valueTot;
  numberOfCharged = valueP;
  if(valueTot < valueP)  { 
    NumberOfExitationWarning("SetNumberOfExcitedParticle"); 
  }
}

inline G4int GXFragment::GetNumberOfHoles()  const
{
  return numberOfHoles;
}

inline G4int GXFragment::GetNumberOfChargedHoles()  const
{
  return numberOfChargedHoles;
}

inline void GXFragment::SetNumberOfHoles(G4int valueTot, G4int valueP)
{
  numberOfHoles = valueTot;
  numberOfChargedHoles = valueP;
  if(valueTot < valueP)  { 
    NumberOfExitationWarning("SetNumberOfHoles"); 
  }
}

inline void GXFragment::SetNumberOfParticles(G4int value)
{
  numberOfParticles = value;
}

inline void GXFragment::SetNumberOfCharged(G4int value)
{
  numberOfCharged = value;
  if(value > numberOfParticles)  { 
    NumberOfExitationWarning("SetNumberOfCharged"); 
  }
}

/*
inline G4int GXFragment::GetNumberOfElectrons() const
{
  return numberOfShellElectrons;
}

inline void GXFragment::SetNumberOfElectrons(G4int value)
{
  numberOfShellElectrons = value;
}
*/

/*
inline G4int GXFragment::GetCreatorModelType() const
{
  return creatorModel;
}

inline void GXFragment::SetCreatorModelType(G4int value)
{
  creatorModel = value;
}
*/

/*
inline G4double GXFragment::GetSpin() const
{
  return spin;
}

inline void GXFragment::SetSpin(G4double value)
{
  spin = value;
}
*/

/*
inline G4int GXFragment::GetFloatingLevelNumber() const
{
  return xLevel;
}

inline void GXFragment::SetFloatingLevelNumber(G4int value)
{
  xLevel = value;
}
*/

/*
inline 
const G4ParticleDefinition* GXFragment::GetParticleDefinition(void) const
{
  return theParticleDefinition;
}

inline void GXFragment::SetParticleDefinition(const G4ParticleDefinition * p)
{
  theParticleDefinition = p;
}
*/

/*
inline G4double GXFragment::GetCreationTime() const 
{
  return theCreationTime;
}

inline void GXFragment::SetCreationTime(G4double time)
{
  theCreationTime = time;
}
*/

/*
inline G4NuclearPolarization* GXFragment::NuclearPolarization()
{
  return thePolarization;
}

inline G4NuclearPolarization* GXFragment::GetNuclearPolarization() const
{
  return thePolarization;
}

inline void GXFragment::SetNuclearPolarization(G4NuclearPolarization* p)
{
  thePolarization = p;
}
*/

#endif
