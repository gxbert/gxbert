#ifndef GXNucleiProperties_h
#define GXNucleiProperties_h 1

#include "globals.hh"
#include "G4ios.hh"

class GXNucleiProperties
{
 // Class Description
 //   GXNucleiProperties is an utility class to provide mass formula of nuclei
 //   (i.e. it has static member function only)

public: 

  // Destructor
  ~GXNucleiProperties() { };

  // Default constructor ()
  GXNucleiProperties() { };


public:  // With Description

  // Give mass of nucleus A,Z
  static G4double GetNuclearMass(const G4double A, const G4double Z);
  static G4double GetNuclearMass(const G4int A, const G4int Z);

  // return 'true' if the nucleus in the stable table 
  // (i.e.in GXNucleiPropertiesTable)
  static bool IsInStableTable(const G4double A, const G4double Z);
  static bool IsInStableTable(const G4int A, const G4int Z);

  // Give binding energy 
  static G4double GetBindingEnergy(const G4int A, const G4int Z);
  static G4double GetBindingEnergy(const G4double A, const G4double Z);

  // Calculate Mass Excess of nucleus A,Z
  static G4double GetMassExcess(const G4int A, const G4int Z);
  static G4double GetMassExcess(const G4double A, const G4double Z);

  //Swich AME table in use
  static void UseOldAMETable( G4bool val = true );

private:
  // hidie methods to enforce using GetNuclearMass
  // Give mass of Atom A,Z
  static G4double GetAtomicMass(const G4double A, const G4double Z);
  
private:
  
  static G4double  AtomicMass(G4double A, G4double Z);
  
  static G4double  NuclearMass(G4double A, G4double Z);
  
  static G4double BindingEnergy(G4double A, G4double Z);
  
  static G4double MassExcess(G4double A, G4double Z);

private: 
  // table of orbit electrons mass - binding energy 
  enum  {MaxZ = 120};
  static G4ThreadLocal G4double electronMass[MaxZ];

private:
  static G4ThreadLocal G4bool   isIntialized;
  static G4ThreadLocal G4double mass_proton;
  static G4ThreadLocal G4double mass_neutron;
  static G4ThreadLocal G4double mass_deuteron;
  static G4ThreadLocal G4double mass_triton;
  static G4ThreadLocal G4double mass_alpha;
  static G4ThreadLocal G4double mass_He3;
  static G4bool use_old_evaluation;
 	
};



#endif
