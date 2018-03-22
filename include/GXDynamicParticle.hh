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
// $Id: GXDynamicParticle.hh 81373 2014-05-27 13:06:32Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//      ---------------- GXDynamicParticle  ----------------
//      first implementation by Makoto Asai, 29 January 1996
//      revised by G.Cosmo, 29 February 1996
//      revised by Hisaya Kurashige, 24 July 1996
//      revised by Hisaya Kurashige, 19 Oct 1996
//      revised by Hisaya Kurashige, 19 Feb 1997
//                ------------------------
//      Add theDynamicCharge and theElectronOccupancy
//                             17 AUg. 1999   H.Kurashige  
//      Add thePreAssignedDecayTime   18 Jan. 2001 H.Kurashige
//      Added  MagneticMoment               Mar. 2007
// ------------------------------------------------------------

#ifndef GXDynamicParticle_h
#define GXDynamicParticle_h 1

#include <cmath>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"

#include "GXParticleDefinition.hh"
#include "GXAllocator.hh"
#include "G4LorentzVector.hh"

//TK for GXBERT
//#include "G4ParticleMomentum.hh"
#include "G4ThreeVector.hh"
typedef G4ThreeVector G4ParticleMomentum;
//  G4ParticleMomentum is "momentum direction" not "momentum vector"
//  The name is miss-leading so you should not use G4ParticleMomentum
//  and you are recommended to use G4ThreeVector instead

#include "GXElectronOccupancy.hh"


//class G4PrimaryParticle;
class  G4VProcess;
class  GXDecayProducts;

class GXDynamicParticle 
{
  // Class Description
  //  The dynamic particle is a class which contains the purely
  //  dynamic aspects of a moving particle. It also has a
  //  pointer to a GXParticleDefinition object, which holds
  //  all the static information.
  //

  public: // With Description
  //- constructors 
     GXDynamicParticle();

     GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
                        const G4ThreeVector& aMomentumDirection,
                        G4double aKineticEnergy);
     GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
                        const G4ThreeVector& aParticleMomentum);
     GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
                        const G4LorentzVector    &aParticleMomentum);
     GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
			G4double aTotalEnergy,
                        const G4ThreeVector &aParticleMomentum);
     GXDynamicParticle(const GXParticleDefinition * aParticleDefinition,
		       const G4ThreeVector& aMomentumDirection,
		       G4double aKineticEnergy,
		       const G4double dynamicalMass);

     GXDynamicParticle(const GXDynamicParticle &right);

  //- destructor
     ~GXDynamicParticle();

  //- operators
     GXDynamicParticle & operator=(const GXDynamicParticle &right);
     G4int operator==(const GXDynamicParticle &right) const;
     G4int operator!=(const GXDynamicParticle &right) const;

  //  new/delete operators are oberloded to use G4Allocator
     inline void *operator new(size_t);
     inline void operator delete(void *aDynamicParticle);

  //- Set/Get methods
 
     const G4ThreeVector& GetMomentumDirection() const;
      //  Returns the normalized direction of the momentum
     void SetMomentumDirection(const G4ThreeVector &aDirection);
      //  Sets the normalized direction of the momentum
     void SetMomentumDirection(G4double px, G4double py, G4double pz);
      //  Sets the normalized direction of the momentum by coordinates

     G4ThreeVector GetMomentum() const;
      //  Returns the current particle momentum vector
     void SetMomentum( const G4ThreeVector &momentum);
      //  set the current particle momentum vector

     G4LorentzVector Get4Momentum() const;
      //  Returns the current particle energy-momentum 4vector
     void Set4Momentum( const G4LorentzVector &momentum);
      //  Set the current particle energy-momentum 4vector


     G4double GetTotalMomentum() const;
      //  Returns the module of the momentum vector
     G4double GetTotalEnergy() const;
      //  Returns the total energy of the particle

     G4double GetKineticEnergy() const;
      //  Returns the kinetic energy of a particle
     void SetKineticEnergy(G4double aEnergy);
      //  Sets the kinetic energy of a particle


     G4double GetProperTime() const;
      //  Returns the current particle proper time
     void SetProperTime( G4double );
      //  Set the current particle Proper Time


     const G4ThreeVector& GetPolarization() const;
     void SetPolarization(G4double polX, G4double polY, G4double polZ);
      //   Set/Get polarization vector       


     G4double GetMass() const;
     void     SetMass(G4double mass);
     // set/get dynamical mass
     // the dynamical mass is set to PDG mass in default


     G4double GetCharge() const;
     void     SetCharge(G4double charge);
     void     SetCharge(G4int    chargeInUnitOfEplus);
     // set/get dynamical charge 
     // the dynamical mass is set to PDG charge in default

     G4double GetSpin() const;
     void     SetSpin(G4double spin);
     void     SetSpin(G4int    spinInUnitOfHalfInteger);
     // set/get dynamical spin
     // the dynamical spin is set to PDG spin in default

     G4double GetMagneticMoment() const;
     void     SetMagneticMoment(G4double magneticMoment);
     // set/get dynamical MagneticMoment  
     // the dynamical mass is set to PDG MagneticMoment in default


     const GXElectronOccupancy* GetElectronOccupancy() const;
     // Get electron occupancy 
     // ElectronOccupancy is valid only if the particle is ion
     G4int  GetTotalOccupancy() const;
     G4int  GetOccupancy(G4int orbit) const;
     void   AddElectron(G4int orbit, G4int number = 1);
     void   RemoveElectron(G4int orbit, G4int number = 1);
 
 
     const GXParticleDefinition* GetParticleDefinition() const;
     void SetDefinition(const GXParticleDefinition * aParticleDefinition);
     //   Set/Get particle definition  
     //  following method of GetDefinition remains 
     //  because of backward compatiblity. It will be removed in future 
     GXParticleDefinition* GetDefinition() const;
 
     
     const GXDecayProducts *GetPreAssignedDecayProducts() const;
     void SetPreAssignedDecayProducts(GXDecayProducts *aDecayProducts);
      //   Set/Get pre-assigned decay channel

     G4double GetPreAssignedDecayProperTime() const;
     void SetPreAssignedDecayProperTime(G4double);
      //   Set/Get pre-assigned proper time when the particle will decay
 
   
  //- print out information
     void DumpInfo(G4int mode= 0) const;
     //    mode 0 : default )(minimum)
     //    mode 1 : 0 + electron occupancy

   protected:
     void      AllocateElectronOccupancy(); 
     G4double  GetElectronMass() const;

   protected:
     G4ThreeVector theMomentumDirection;
      //  The normalized momentum vector

     const GXParticleDefinition *theParticleDefinition;
      //  Contains the static information of this particle.

     G4ThreeVector thePolarization;

     G4double theKineticEnergy;

     G4double theProperTime;

     G4double theDynamicalMass;

     G4double theDynamicalCharge;

     G4double theDynamicalSpin;

     G4double theDynamicalMagneticMoment;

     GXElectronOccupancy* theElectronOccupancy;          
  
     GXDecayProducts *thePreAssignedDecayProducts;

     G4double thePreAssignedDecayTime;

 protected:
   G4int verboseLevel;
 
 public:  // With Description
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // Set/Get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

 protected:
   //G4PrimaryParticle* primaryParticle;
   // This void pointer is used by G4EventManager to maintain the
   // link between pre-assigned decay products and corresponding
   // primary particle.

 public:
   //void SetPrimaryParticle(G4PrimaryParticle* p);
   void SetPDGcode(G4int c);

 public: // With Description
   //G4PrimaryParticle* GetPrimaryParticle() const;
   // Return the pointer to the corresponding G4PrimaryParticle object
   // if this particle is a primary particle OR is defined as a pre-assigned
   // decay product. Otherwise return null.

   G4int GetPDGcode() const;
   // Return the PDG code of this particle. If the particle is known to Geant4
   // its PDG code defined in GXParticleDefinition is returned. If it is unknown
   // (i.e. PDG code in GXParticleDefinition is 0), PDG code defined in the
   // corresponding primary particle or pre-assigned decay product will be
   // returned if available. Otherwise (e.g. for geantino) returns 0.

 protected:
   G4int thePDGcode;
};

#include "GXDynamicParticle.icc"

#endif
