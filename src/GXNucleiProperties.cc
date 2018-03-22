#include "GXNucleiProperties.hh"

#include "GXNucleiPropertiesTableAME12.hh"
#include "GXNucleiPropertiesTheoreticalTable.hh"
#include "GXParticleTable.hh"

#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"

G4ThreadLocal G4double GXNucleiProperties::mass_proton = -1.;
G4ThreadLocal G4double GXNucleiProperties::mass_neutron = -1.;
G4ThreadLocal G4double GXNucleiProperties::mass_deuteron = -1.;
G4ThreadLocal G4double GXNucleiProperties::mass_triton = -1.;
G4ThreadLocal G4double GXNucleiProperties::mass_alpha = -1.;
G4ThreadLocal G4double GXNucleiProperties::mass_He3 = -1.;
#ifndef GXNucleiProperties_USE_OLD_AME_TABLE
G4bool GXNucleiProperties::use_old_evaluation = false;
#else
G4bool GXNucleiProperties::use_old_evaluation = true;
#endif

G4double GXNucleiProperties::GetNuclearMass(const G4double A, const G4double Z)
{
  G4double mass=0.0;

  if (std::fabs(A - G4int(A)) > 1.e-10) {
    mass = NuclearMass(A,Z);
 
  } else {
    // use mass table
    G4int iZ = G4int(Z);
    G4int iA = G4int(A);
    mass =GetNuclearMass(iA,iZ);
  }
  
   return mass;
}


G4double GXNucleiProperties::GetNuclearMass(const G4int A, const G4int Z)
{
  if (mass_proton  <= 0.0 ) {
    const GXParticleDefinition * nucleus = 0;
    nucleus = GXParticleTable::GetParticleTable()->FindParticle("proton"); // proton 
    if (nucleus!=0) mass_proton = nucleus->GetPDGMass();
    nucleus = GXParticleTable::GetParticleTable()->FindParticle("neutron"); // neutron 
    if (nucleus!=0) mass_neutron = nucleus->GetPDGMass();
    nucleus = GXParticleTable::GetParticleTable()->FindParticle("deuteron"); // deuteron 
    if (nucleus!=0) mass_deuteron = nucleus->GetPDGMass();
    nucleus = GXParticleTable::GetParticleTable()->FindParticle("triton"); // triton 
    if (nucleus!=0) mass_triton = nucleus->GetPDGMass();
    nucleus = GXParticleTable::GetParticleTable()->FindParticle("alpha"); // alpha 
    if (nucleus!=0) mass_alpha = nucleus->GetPDGMass();
    nucleus = GXParticleTable::GetParticleTable()->FindParticle("He3"); // He3 
    if (nucleus!=0) mass_He3 = nucleus->GetPDGMass();

  }

  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (GXParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "GXNucleiProperties::GetNuclearMass: Wrong values for A = " << A 
	     << " and Z = " << Z << G4endl;
    }
#endif    
    return 0.0;
  }
  
  G4double mass= -1.;
  if ( (Z<=2) ) {
    // light nuclei
    if ( (Z==1)&&(A==1) ) {
      mass = mass_proton;
    } else if ( (Z==0)&&(A==1) ) {
      mass = mass_neutron;
    } else if ( (Z==1)&&(A==2) ) {
      mass = mass_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      mass = mass_triton;
    } else if ( (Z==2)&&(A==4) ) {
      mass = mass_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      mass = mass_He3;
    }
  }
  
  if (mass < 0.) {
    G4bool inAMETable = false;
    //if ( ! use_old_evaluation ) {
       inAMETable = GXNucleiPropertiesTableAME12::IsInTable(Z,A);
    //} else {
    //   inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    //}
    if ( inAMETable ) {
      // AME table
      //if ( ! use_old_evaluation ) {
         mass = GXNucleiPropertiesTableAME12::GetNuclearMass(Z,A);
      //} else {
      //   mass = G4NucleiPropertiesTableAME03::GetNuclearMass(Z,A);
      //}
    } else if (GXNucleiPropertiesTheoreticalTable::IsInTable(Z,A)){
      // Theoretical table
      mass = GXNucleiPropertiesTheoreticalTable::GetNuclearMass(Z,A);
    } else {
      mass = NuclearMass(G4double(A),G4double(Z));
    }
  }

  if (mass < 0.) mass = 0.0;
  return mass;
}

G4bool GXNucleiProperties::IsInStableTable(const G4double A, const G4double Z)
{
  G4int iA = G4int(A);
  G4int iZ = G4int(Z);
  return IsInStableTable(iA, iZ);
}

G4bool GXNucleiProperties::IsInStableTable(const G4int A, const int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (GXParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "GXNucleiProperties::IsInStableTable: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;	
    }
#endif 
    return false;
  } 

  //if ( ! use_old_evaluation ) {
    return GXNucleiPropertiesTableAME12::IsInTable(Z,A);
  //} else {
  //  return GXNucleiPropertiesTableAME03::IsInTable(Z,A);
  //}
}

G4double GXNucleiProperties::GetMassExcess(const G4double A, const G4double Z)
{
  G4int iA = G4int(A);
  G4int iZ = G4int(Z);
  return GetMassExcess(iA,iZ);
}

G4double GXNucleiProperties::GetMassExcess(const G4int A, const G4int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (GXParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "GXNucleiProperties::GetMassExccess: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;
    }
#endif    
    return 0.0;
    
  } else {

    G4bool inAMETable = false;
    //if ( ! use_old_evaluation ) {
       inAMETable = GXNucleiPropertiesTableAME12::IsInTable(Z,A);
    //} else {
    //   inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    //}
    if (inAMETable){
      // AME table
      //if ( ! use_old_evaluation ) {
         return GXNucleiPropertiesTableAME12::GetMassExcess(Z,A);
      //} else {
      //   return G4NucleiPropertiesTableAME03::GetMassExcess(Z,A);
      //}
    } else if (GXNucleiPropertiesTheoreticalTable::IsInTable(Z,A)){
      return GXNucleiPropertiesTheoreticalTable::GetMassExcess(Z,A);
    } else {
      return MassExcess(A,Z);
    }
  }

}


G4double GXNucleiProperties::GetAtomicMass(const G4double A, const G4double Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (GXParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "GXNucleiProperties::GetAtomicMass: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;	
    }
#endif 
    return 0.0;

  } else if (std::fabs(A - G4int(A)) > 1.e-10) {
    return AtomicMass(A,Z);

  } else {
    G4int iA = G4int(A);
    G4int iZ = G4int(Z);
    G4bool inAMETable = false;
    //if ( ! use_old_evaluation ) {
       inAMETable = GXNucleiPropertiesTableAME12::IsInTable(Z,A);
    //} else {
    //   inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    //}
    if (inAMETable) {
      //if ( ! use_old_evaluation ) {
        return GXNucleiPropertiesTableAME12::GetAtomicMass(Z,A);
      //} else {
      //  return G4NucleiPropertiesTableAME03::GetAtomicMass(Z,A);
      //}
    } else if (GXNucleiPropertiesTheoreticalTable::IsInTable(iZ,iA)){
      return GXNucleiPropertiesTheoreticalTable::GetAtomicMass(iZ,iA);
    } else {
      return AtomicMass(A,Z);
    }
  }
}

G4double GXNucleiProperties::GetBindingEnergy(const G4double A, const G4double Z)
{
  G4int iA = G4int(A);
  G4int iZ = G4int(Z);
  return GetBindingEnergy(iA,iZ);
}

G4double GXNucleiProperties::GetBindingEnergy(const G4int A, const G4int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (GXParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "GXNucleiProperties::GetMassExccess: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;
    }
#endif
    return 0.0;

  } else {
    G4bool inAMETable = false;
    //if ( ! use_old_evaluation ) {
       inAMETable = GXNucleiPropertiesTableAME12::IsInTable(Z,A);
    //} else {
    //   inAMETable = G4NucleiPropertiesTableAME03::IsInTable(Z,A);
    //}
    if (inAMETable) {
      //if ( ! use_old_evaluation ) {
        return GXNucleiPropertiesTableAME12::GetBindingEnergy(Z,A);
      //} else {
      //  return G4NucleiPropertiesTableAME03::GetBindingEnergy(Z,A);
      //}
    } else if (GXNucleiPropertiesTheoreticalTable::IsInTable(Z,A)) {
      return GXNucleiPropertiesTheoreticalTable::GetBindingEnergy(Z,A);
    }else {
      return BindingEnergy(A,Z);
    }

  }
}


G4double GXNucleiProperties::MassExcess(G4double A, G4double Z) 
{
  return GetAtomicMass(A,Z) - A*amu_c2;
}

G4double  GXNucleiProperties::AtomicMass(G4double A, G4double Z)
{
  //const G4double hydrogen_mass_excess;
  //const G4double neutron_mass_excess;  
  G4double hydrogen_mass_excess;
  G4double neutron_mass_excess;  
  //if ( ! use_old_evaluation ) {
    hydrogen_mass_excess = GXNucleiPropertiesTableAME12::GetMassExcess(1,1);
    neutron_mass_excess = GXNucleiPropertiesTableAME12::GetMassExcess(0,1);
  //} else {
  //  hydrogen_mass_excess = G4NucleiPropertiesTableAME03::GetMassExcess(1,1);
  //  neutron_mass_excess = G4NucleiPropertiesTableAME03::GetMassExcess(0,1);
  //}

  G4double mass =
      (A-Z)*neutron_mass_excess + Z*hydrogen_mass_excess - BindingEnergy(A,Z) + A*amu_c2;

  return mass;
}

G4double  GXNucleiProperties::NuclearMass(G4double A, G4double Z)
{
  if (A < 1 || Z < 0 || Z > A) {
#ifdef G4VERBOSE
    if (GXParticleTable::GetParticleTable()->GetVerboseLevel()>0) {
      G4cout << "GXNucleiProperties::NuclearMass: Wrong values for A = " 
	     << A << " and Z = " << Z << G4endl;
    }
#endif 
    return 0.0;
  }

  G4double mass = AtomicMass(A,Z);
  // atomic mass is converted to nuclear mass according formula in  AME03 and 12
  mass -= Z*electron_mass_c2;
  mass += ( 14.4381*std::pow ( Z , 2.39 ) + 1.55468*1e-6*std::pow ( Z , 5.35 ) )*eV;      

  return mass;
}

G4double  GXNucleiProperties::BindingEnergy(G4double A, G4double Z)
{ 
  //
  // Weitzsaecker's Mass formula
  //
  G4int Npairing = G4int(A-Z)%2;                  // pairing
  G4int Zpairing = G4int(Z)%2;
  G4double binding =
      - 15.67*A                           // nuclear volume
      + 17.23*std::pow(A,2./3.)                // surface energy
      + 93.15*((A/2.-Z)*(A/2.-Z))/A       // asymmetry
      + 0.6984523*Z*Z*std::pow(A,-1./3.);      // coulomb
  if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / std::sqrt(A);  // pairing

  return -binding*MeV;
}

void GXNucleiProperties::UseOldAMETable( G4bool val )
{
   use_old_evaluation = val;
}
