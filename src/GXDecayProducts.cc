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
// $Id: GXDecayProducts.cc 69015 2013-04-15 09:46:48Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      10 July 1996 H.Kurashige
//      21 Oct  1996 H.Kurashige
//      12 Dec 1997 H.Kurashige
//      4 Apr. 2012 H.Kurashige use std::vector
// ------------------------------------------------------------

#include "G4ios.hh"
#include "globals.hh"
#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"
#include "GXDecayProducts.hh"

#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


GXDecayProducts::GXDecayProducts()
                :numberOfProducts(0),theParentParticle(0)
{ 
   theProductVector = new G4DecayProductVector();
}

GXDecayProducts::GXDecayProducts(const GXDynamicParticle &aParticle)
  :numberOfProducts(0),theParentParticle(0)
{
  theParentParticle = new GXDynamicParticle(aParticle);
  theProductVector = new G4DecayProductVector();
}

GXDecayProducts::GXDecayProducts(const GXDecayProducts &right) 
                :numberOfProducts(0)
{
  theProductVector = new G4DecayProductVector();

  // copy parent (Deep Copy)
  theParentParticle = new GXDynamicParticle(*right.theParentParticle);

  //copy daughters (Deep Copy)
  for (G4int index=0; index < right.numberOfProducts; index++) {
    GXDynamicParticle* daughter = right.theProductVector->at(index);
    GXDynamicParticle* pDaughter =  new GXDynamicParticle(*daughter);

    G4double properTime = daughter->GetPreAssignedDecayProperTime();
    if(properTime>0.0)pDaughter->SetPreAssignedDecayProperTime(properTime); 

    //const GXDecayProducts* pPreAssigned = daughter->GetPreAssignedDecayProducts();
    //if (pPreAssigned) {
    //  GXDecayProducts* pPA = new GXDecayProducts(*pPreAssigned);
    //  pDaughter->SetPreAssignedDecayProducts(pPA);
    //}
    //TK Modefied above 5 lines
    if ( daughter->GetPreAssignedDecayProducts() ) {
       G4cout << "GXDecayProducts does not support PreAssignedDecayProducts" << G4endl;
       abort();
    }

    theProductVector->push_back( pDaughter );
  }
  numberOfProducts = right.numberOfProducts;
}

GXDecayProducts & GXDecayProducts::operator=(const GXDecayProducts &right)
{
  G4int index;

  if (this != &right)
  { 
    // recreate parent
    if (theParentParticle != 0) delete theParentParticle;
    theParentParticle = new GXDynamicParticle(*right.theParentParticle);

    // delete GXDynamicParticle objects
    for (index=0; index < numberOfProducts; index++) {
      delete theProductVector->at(index);
    }
    theProductVector->clear();

    //copy daughters (Deep Copy)
    for (index=0; index < right.numberOfProducts; index++) {
      GXDynamicParticle* daughter = right.theProductVector->at(index);
      GXDynamicParticle* pDaughter =  new GXDynamicParticle(*daughter);

      G4double properTime = daughter->GetPreAssignedDecayProperTime();
      if(properTime>0.0) pDaughter->SetPreAssignedDecayProperTime(properTime); 
      
//    const GXDecayProducts* pPreAssigned = daughter->GetPreAssignedDecayProducts();
//    if (pPreAssigned) {
//  	GXDecayProducts* pPA = new GXDecayProducts(*pPreAssigned);
//	pDaughter->SetPreAssignedDecayProducts(pPA);
//    }
      //TK Modefied above 5 lines
      if ( daughter->GetPreAssignedDecayProducts() ) {
         G4cout << "GXDecayProducts does not support PreAssignedDecayProducts" << G4endl;
         abort();
      }
      theProductVector->push_back( pDaughter );
    } 
    numberOfProducts = right.numberOfProducts;
    
  }
  return *this;
}

GXDecayProducts::~GXDecayProducts()
{
  //delete parent
  if (theParentParticle != 0) delete theParentParticle;
  
  // delete GXDynamicParticle object
  for (G4int index=0; index < numberOfProducts; index++) {
      delete theProductVector->at(index);
  }
  theProductVector->clear();
  numberOfProducts = 0;    
  delete theProductVector;
}

GXDynamicParticle* GXDecayProducts::PopProducts()
{
   if ( numberOfProducts >0 ) {
     numberOfProducts -= 1;   
     GXDynamicParticle* part = theProductVector->back();
     theProductVector->pop_back();
     return part;
   } else {
     return 0;
   }
}

G4int GXDecayProducts::PushProducts(GXDynamicParticle *aParticle)
{
  theProductVector->push_back(aParticle);
  numberOfProducts += 1; 
  return numberOfProducts;
}

GXDynamicParticle* GXDecayProducts::operator[](G4int anIndex) const
{
   if ((numberOfProducts > anIndex) && (anIndex >=0) ) {
     return  theProductVector->at(anIndex);
   } else {
     return 0;
   }
}

void  GXDecayProducts::SetParentParticle(const GXDynamicParticle &aParticle)
{
  if (theParentParticle != 0) delete theParentParticle;
  theParentParticle = new GXDynamicParticle(aParticle);
}


void GXDecayProducts::Boost(G4double totalEnergy, const G4ThreeVector &momentumDirection)
{
  // calcurate new beta
  G4double   mass = theParentParticle->GetMass();
  G4double   totalMomentum(0);
  if (totalEnergy > mass ) totalMomentum  = std::sqrt( (totalEnergy - mass)*(totalEnergy + mass) );
  G4double   betax = momentumDirection.x()*totalMomentum/totalEnergy;  
  G4double   betay = momentumDirection.y()*totalMomentum/totalEnergy;  
  G4double   betaz = momentumDirection.z()*totalMomentum/totalEnergy;  
  this->Boost(betax, betay, betaz);
}

void GXDecayProducts::Boost(G4double newbetax, G4double newbetay, G4double newbetaz)
{ 
  G4double   mass = theParentParticle->GetMass();
  G4double   energy  = theParentParticle->GetTotalEnergy();
  G4double   momentum  = 0.0;

  G4ThreeVector direction(0.0,0.0,1.0);    
  G4LorentzVector p4;

  if (energy - mass > DBL_MIN) {
    // calcurate  beta of initial state
    momentum  = theParentParticle->GetTotalMomentum();
    direction = theParentParticle->GetMomentumDirection();
    G4double betax = -1.0*direction.x()*momentum/energy;  
    G4double betay = -1.0*direction.y()*momentum/energy;  
    G4double betaz = -1.0*direction.z()*momentum/energy;  
    
    for (G4int index=0; index < numberOfProducts; index++) {
       // make G4LorentzVector for secondaries
       p4 = (theProductVector->at(index))->Get4Momentum();

       // boost secondaries to theParentParticle's rest frame 
       p4.boost(betax, betay, betaz);

       // boost secondaries to  new frame 
       p4.boost(newbetax, newbetay, newbetaz);

       // change energy/momentum
       (theProductVector->at(index))->Set4Momentum(p4);
     }
   } else {
     for (G4int index=0; index < numberOfProducts; index++) {
       // make G4LorentzVector for secondaries
       p4 = (theProductVector->at(index))->Get4Momentum();

       // boost secondaries to  new frame 
       p4.boost(newbetax, newbetay, newbetaz);

       // change energy/momentum
       (theProductVector->at(index))->Set4Momentum(p4);
      }
   }
   // make G4LorentzVector for parent in its rest frame
   mass = theParentParticle->GetMass();
   G4LorentzVector parent4( 0.0, 0.0, 0.0, mass);

   // boost parent to new frame 
   parent4.boost(newbetax, newbetay, newbetaz);

   // change energy/momentum
   theParentParticle->Set4Momentum(parent4);
}

G4bool GXDecayProducts::IsChecked() const
{
  G4bool returnValue = true;
  // check parent 
  //   energy/momentum
  G4double   parent_energy  = theParentParticle->GetTotalEnergy();
  G4ThreeVector direction = theParentParticle->GetMomentumDirection();
  G4ThreeVector parent_momentum = direction*(theParentParticle->GetTotalMomentum());
  // check momentum dirction is a unit vector
  if ( (parent_momentum.mag() >0.0) && (std::fabs(direction.mag()-1.0) >1.0e-6 ) ) {
#ifdef G4VERBOSE
    G4cout << "GXDecayProducts::IsChecked()::  "
           << " Momentum Direction Vector of Parent is not normalized "
           << "  (=" << direction.mag() << ")" << G4endl;
#endif
    returnValue = false;
    parent_momentum = parent_momentum * (1./direction.mag());
  }

  //daughters
  G4double   mass, energy;
  G4ThreeVector momentum;
  G4double   total_energy = parent_energy;
  G4ThreeVector total_momentum =  parent_momentum;
  for (G4int index=0; index < numberOfProducts; index++) 
  {
    GXDynamicParticle* part = theProductVector->at(index);
    mass = part->GetMass();
    energy  = part->GetTotalEnergy();
    direction = part->GetMomentumDirection();
    momentum = direction*(part->GetTotalMomentum());
    // check momentum dirction is a unit vector
    if ( (momentum.mag()>0.0) && (std::fabs(direction.mag()-1.0) > 1.0e-6)) {
#ifdef G4VERBOSE
      G4cout <<  "GXDecayProducts::IsChecked()::  "
             << " Momentum Direction Vector of Daughter [" << index
             << "]  is not normalized (=" << direction.mag() << ")" << G4endl;
#endif
      returnValue = false;
      momentum = momentum * (1./direction.mag());
    }
    // whether daughter stops or not
    if (energy - mass < DBL_MIN ) {
#ifdef G4VERBOSE
      G4cout <<  "GXDecayProducts::IsChecked()::  "
             << "  Daughter [" << index << "] has no kinetic energy "<< G4endl;
#endif
      returnValue = false;
    }
    total_energy -= energy; 
    total_momentum -= momentum;
  }
  // check energy/momentum conservation
  if ( (std::fabs(total_energy) >1.0e-9*MeV) || (total_momentum.mag() >1.0e-9*MeV ) ){ 
#ifdef G4VERBOSE
    G4cout <<  "GXDecayProducts::IsChecked()::  "
           << " Energy/Momentum is not conserved   "<< G4endl;
    G4cout << " difference between parent energy and sum of dughters' energy : " 
	   << total_energy /MeV << "[MeV]  " << G4endl; 
    G4cout << " difference between parent momentum and sum of dughters' momentum : " 
	   << " x:" << total_momentum.getX()/MeV 
	   << " y:" << total_momentum.getY()/MeV  
	   << " z:" << total_momentum.getZ()/MeV  
	   << G4endl;
#endif
    returnValue = false;
  }
  return returnValue;
}

void GXDecayProducts::DumpInfo() const
{
   G4cout << " ----- List of DecayProducts  -----" << G4endl;
   G4cout << " ------ Parent Particle ----------" << G4endl;
   if (theParentParticle != 0) theParentParticle->DumpInfo();
   G4cout << " ------ Daughter Particles  ------" << G4endl;  
   for (G4int index=0; index < numberOfProducts; index++) 
   {
      G4cout << " ----------" << index+1 << " -------------" << G4endl;  
      (theProductVector->at(index))-> DumpInfo();
   }
   G4cout << " ----- End List of DecayProducts  -----" << G4endl;
   G4cout << G4endl;
} 
