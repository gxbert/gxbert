//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: GXDynamicParticle.cc 98352 2016-07-08 08:21:00Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//      ---------------- GXDynamicParticle  ----------------
//      first implementation by Makoto Asai, 29 January 1996
//      revised by G.Cosmo, 29 February 1996
//      revised by H.Kurashige 06 May 1996
//      revised by Hisaya Kurashige, 27 July 1996
//         modify thePreAssignedDecayProducts
//         add   void SetMomentum(G4ThreeVector &momentum)
//         add   void Set4Momentum(G4LorentzVector &momentum)
//         add   GXDynamicParticle(GXParticleDefinition * aParticleDefinition,
//                                 G4LorentzVector &p4vector)
//      revised by Hisaya Kurashige, 19 Oct 1996
//         add    theKillProcess
//         add    ProperTime
//      revised by Hisaya Kurashige, 26 Mar 1997
//         modify destructor
//      revised by Hisaya Kurashige, 05 June 1997
//         modify DumpInfo()
//      revised by Hisaya Kurashige, 5  June 1998
//         remove    theKillProcess
//      revised by Hisaya Kurashige, 5  Mar 2001
//         fixed  SetDefinition()
//      revised by V.Ivanchenko,    12 June 2003
//         fixed problem of massless particles
//      revised by V.Ivanchenko,    18 June 2003
//         take into account the case of virtual photons
//      revised by M.Kelsey         12 May 2010
//	   ensure that all constructors initialize all data members
//--------------------------------------------------------------

#include "GXDynamicParticle.hh"
#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"
#include "GXDecayProducts.hh"
#include "G4LorentzVector.hh"
#include "GXParticleDefinition.hh"
#include "GXParticleTable.hh"
#include "GXIonTable.hh"
//#include "G4PrimaryParticle.hh"

G4ThreadLocal GXAllocator<GXDynamicParticle> *pGXDynamicParticleAllocator = 0;

static const G4double EnergyMomentumRelationAllowance = 1.0e-2*keV;

////////////////////
GXDynamicParticle::GXDynamicParticle():
		   theMomentumDirection(0.0,0.0,1.0),
		   theParticleDefinition(0),
		   theKineticEnergy(0.0),
 		   theProperTime(0.0),
		   theDynamicalMass(0.0),
		   theDynamicalCharge(0.0),
		   theDynamicalSpin(0.0),
		   theDynamicalMagneticMoment(0.0),
		   theElectronOccupancy(0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
		   //primaryParticle(0),
                   thePDGcode(0)
{
}

////////////////////
// -- constructors ----
////////////////////
GXDynamicParticle::GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
				     const G4ThreeVector& aMomentumDirection,
				     G4double aKineticEnergy):
		   theMomentumDirection(aMomentumDirection),
		   theParticleDefinition(aParticleDefinition),
		   theKineticEnergy(aKineticEnergy),
 		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
		   theElectronOccupancy(0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
		   //primaryParticle(0),
                   thePDGcode(0)
{
}

////////////////////
GXDynamicParticle::GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
				     const G4ThreeVector& aMomentumDirection,
				     G4double aKineticEnergy,
                                     const G4double dynamicalMass):
		   theMomentumDirection(aMomentumDirection),
		   theParticleDefinition(aParticleDefinition),
		   theKineticEnergy(aKineticEnergy),
 		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
		   theElectronOccupancy(0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
		   //primaryParticle(0),
                   thePDGcode(0)
{
  if (std::abs(theDynamicalMass-dynamicalMass)> EnergyMomentumRelationAllowance) {
    if (dynamicalMass>EnergyMomentumRelationAllowance) theDynamicalMass= dynamicalMass;
    else  theDynamicalMass= 0.0;
  } 
}

////////////////////
GXDynamicParticle::GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
                                     const G4ThreeVector& aParticleMomentum):
		   theParticleDefinition(aParticleDefinition),
		   theKineticEnergy(0.0),
       		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
		   theElectronOccupancy(0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
		   //primaryParticle(0),
                   thePDGcode(0)
{
  SetMomentum(aParticleMomentum);  // 3-dim momentum is given
}


////////////////////
GXDynamicParticle::GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
				     const G4LorentzVector   &aParticleMomentum):
		   theParticleDefinition(aParticleDefinition),
		   theKineticEnergy(0.0),
 		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
		   theElectronOccupancy(0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
		   //primaryParticle(0),
                   thePDGcode(0)
{
  Set4Momentum(aParticleMomentum);  // 4-momentum vector (Lorentz vector) is given
}

GXDynamicParticle::GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
                                     G4double totalEnergy,
				     const G4ThreeVector &aParticleMomentum):
                   theParticleDefinition(aParticleDefinition),
		   theKineticEnergy(0.0),
                   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
		   theElectronOccupancy(0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
		   //primaryParticle(0),
                   thePDGcode(0)
{
  // total energy and 3-dim momentum are given
  G4double pModule2 = aParticleMomentum.mag2();
  if (pModule2>0.0) {
    G4double mass2 = totalEnergy*totalEnergy - pModule2;
    G4double PDGmass2 = (aParticleDefinition->GetPDGMass())*(aParticleDefinition->GetPDGMass());
    SetMomentumDirection(aParticleMomentum.unit());
    if (mass2 < EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance) {
      theDynamicalMass = 0.;
      SetKineticEnergy(totalEnergy);
    } else {
      if (std::abs(PDGmass2-mass2)>EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance){
      theDynamicalMass = std::sqrt(mass2);
      SetKineticEnergy(totalEnergy-theDynamicalMass);
      } else {
	SetKineticEnergy(totalEnergy-theDynamicalMass);
      }
    }
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

////////////////////
GXDynamicParticle::GXDynamicParticle(const GXDynamicParticle &right):
  theMomentumDirection(right.theMomentumDirection),
  theParticleDefinition(right.theParticleDefinition),
  thePolarization(right.thePolarization),
  theKineticEnergy(right.theKineticEnergy),
  theProperTime(0.0),
  theDynamicalMass(right.theDynamicalMass),
  theDynamicalCharge(right.theDynamicalCharge),
  theDynamicalSpin(right.theDynamicalSpin),
  theDynamicalMagneticMoment(right.theDynamicalMagneticMoment),
  theElectronOccupancy(0),
  thePreAssignedDecayProducts(0),	// Do not copy preassignedDecayProducts
  thePreAssignedDecayTime(-1.0),
  verboseLevel(right.verboseLevel),
  //primaryParticle(right.primaryParticle),
  thePDGcode(right.thePDGcode)
{
  if (right.theElectronOccupancy != 0) {
      theElectronOccupancy =
	new GXElectronOccupancy(*right.theElectronOccupancy);
  }
}

////////////////////
// -- destructor ----
////////////////////
GXDynamicParticle::~GXDynamicParticle()
{
  //  delete thePreAssignedDecayProducts
  if (thePreAssignedDecayProducts != 0) delete thePreAssignedDecayProducts;
  thePreAssignedDecayProducts = 0;

  if (theElectronOccupancy != 0) delete theElectronOccupancy;
  theElectronOccupancy =0;
}


////////////////////
// -- operators ----
////////////////////
GXDynamicParticle & GXDynamicParticle::operator=(const GXDynamicParticle &right)
{
  if (this != &right) {
    theMomentumDirection = right.theMomentumDirection;
    theParticleDefinition = right.theParticleDefinition;
    thePolarization = right.thePolarization;
    theKineticEnergy = right.theKineticEnergy;
    theProperTime = right.theProperTime;

    theDynamicalMass = right.theDynamicalMass;
    theDynamicalCharge = right.theDynamicalCharge;
    theDynamicalSpin = right.theDynamicalSpin;
    theDynamicalMagneticMoment = right.theDynamicalMagneticMoment;

    if (theElectronOccupancy != 0) delete theElectronOccupancy;
    if (right.theElectronOccupancy != 0){
      theElectronOccupancy =
             new GXElectronOccupancy(*right.theElectronOccupancy);
    } else {
      theElectronOccupancy = 0;
    }

    // thePreAssignedDecayProducts must not be copied.
    thePreAssignedDecayProducts = 0;
    thePreAssignedDecayTime = -1.0;

    verboseLevel = right.verboseLevel;

    // Primary particle information must be preserved
    //*** primaryParticle = right.primaryParticle;

    thePDGcode = right.thePDGcode;
  }
  return *this;
}

////////////////////
void GXDynamicParticle::SetDefinition(const GXParticleDefinition * aParticleDefinition)
{
  // remove preassigned decay
  if (thePreAssignedDecayProducts != 0) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " GXDynamicParticle::SetDefinition()::"
             << "!!! Pre-assigned decay products is attached !!!! " << G4endl;
      G4cout << "!!! New Definition is " << aParticleDefinition->GetParticleName() 
	     << " !!! " << G4endl;
      G4cout << "!!! Pre-assigned decay products will be deleted !!!! " << G4endl;
    }
#endif
    delete thePreAssignedDecayProducts;
  }
  thePreAssignedDecayProducts = 0;

  theParticleDefinition = aParticleDefinition;

  // set Dynamic mass/chrge
  theDynamicalMass = theParticleDefinition->GetPDGMass();
  theDynamicalCharge = theParticleDefinition->GetPDGCharge();
  theDynamicalSpin = theParticleDefinition->GetPDGSpin();
  theDynamicalMagneticMoment = theParticleDefinition->GetPDGMagneticMoment();

  // Set electron orbits
  if (theElectronOccupancy != 0) delete theElectronOccupancy;
  theElectronOccupancy =0;
  //AllocateElectronOccupancy();

}

////////////////////
G4int GXDynamicParticle::operator==(const GXDynamicParticle &right) const
{
  return (this == (GXDynamicParticle *) &right);
}

////////////////////
G4int GXDynamicParticle::operator!=(const GXDynamicParticle &right) const
{
  return (this != (GXDynamicParticle *) &right);
}



////////////////////
// -- AllocateElectronOccupancy --
////////////////////
void  GXDynamicParticle::AllocateElectronOccupancy()
{
  const GXParticleDefinition* particle = GetDefinition();

  if (GXIonTable::IsIon(particle)) {
    // Only ions can have ElectronOccupancy
    theElectronOccupancy = new GXElectronOccupancy();

  } else {
    theElectronOccupancy = 0;

  }
}

////////////////////
// -- methods for setting Energy/Momentum  --
////////////////////
void GXDynamicParticle::SetMomentum(const G4ThreeVector &momentum)
{
  G4double pModule2 = momentum.mag2();
  if (pModule2>0.0) {
    G4double mass = theDynamicalMass;
    SetMomentumDirection(momentum.unit());
    SetKineticEnergy(std::sqrt(pModule2 + mass*mass)-mass);
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

////////////////////
void GXDynamicParticle::Set4Momentum(const G4LorentzVector &momentum )
{
  G4double pModule2 = momentum.vect().mag2();
  if (pModule2>0.0) {
    SetMomentumDirection(momentum.vect().unit());
    G4double totalenergy = momentum.t();
    G4double mass2 = totalenergy*totalenergy - pModule2;
    G4double PDGmass2 = (theParticleDefinition->GetPDGMass())*(theParticleDefinition->GetPDGMass());
    if (mass2 < EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance) {
      theDynamicalMass = 0.;
    } else if (std::abs(PDGmass2-mass2)>EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance){
      theDynamicalMass = std::sqrt(mass2);
    }
    SetKineticEnergy(totalenergy-theDynamicalMass);
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}


////////////////////
//  --- Dump Information --
////////////////////
#ifdef G4VERBOSE
void GXDynamicParticle::DumpInfo(G4int mode) const
{
  if (theParticleDefinition == 0) {
    G4cout << " GXDynamicParticle::DumpInfo():: !!!Particle type not defined !!!! " << G4endl;
  } else {
    G4cout << " Particle type - " << theParticleDefinition->GetParticleName() << G4endl
         << "   mass:        " << GetMass()/GeV <<  "[GeV]" <<G4endl
         << "   charge:      " << GetCharge()/eplus <<  "[e]" <<G4endl
         << "   Direction x: " << GetMomentumDirection().x() << ", y: "
	 << GetMomentumDirection().y() << ", z: "
                             << GetMomentumDirection().z() << G4endl
         << "   Total Momentum = " << GetTotalMomentum() /GeV << "[GeV]" << G4endl
         << "   Momentum: "    << GetMomentum().x() /GeV << "[GeV]" << ", y: "
                               << GetMomentum().y() /GeV << "[GeV]" << ", z: "
                               << GetMomentum().z() /GeV << "[GeV]" << G4endl
         << "   Total Energy   = " << GetTotalEnergy()/GeV << "[GeV]"  << G4endl
         << "   Kinetic Energy = " << GetKineticEnergy() /GeV << "[GeV]" << G4endl
	 << " MagneticMoment  [MeV/T]: " << GetMagneticMoment()/MeV*tesla << G4endl
         << "   ProperTime     = " << GetProperTime() /ns <<  "[ns]" << G4endl;

    if (mode>0) {
      if( theElectronOccupancy != 0) {
	theElectronOccupancy->DumpInfo();
      }
    }
  }
}
#else
void GXDynamicParticle::DumpInfo(G4int) const
{
  return;
}
#endif
///////////////////////
//
void GXDynamicParticle::SetMass(G4double newMass)
{
  if(std::abs(newMass-theParticleDefinition->GetPDGMass())>EnergyMomentumRelationAllowance*0.1) {
    theDynamicalMass = newMass;
  }
}
////////////////////////
G4double  GXDynamicParticle::GetElectronMass() const
{
  static G4ThreadLocal G4double electronMass = 0.0;

  // check if electron exits and get the mass
  if (electronMass<=0.0) {
    GXParticleDefinition* electron = GXParticleTable::GetParticleTable()->FindParticle("e-");
    if (electron == 0) {
      G4Exception("GXDynamicParticle::GetElectronMass()","PART021",
 		  FatalException,"GXDynamicParticle: GXElectron is not defined !!");
    }
    electronMass = electron->GetPDGMass();
  }

  return electronMass;
}
