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
// $Id: GXDecayTable.cc 105720 2017-08-16 12:38:10Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
// ----------------------------------------
//      implementation for STL          14 Feb. 2000 H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "GXDecayTable.hh"
#include "Randomize.hh"

GXDecayTable::GXDecayTable()
  :parent(nullptr), channels(nullptr)
{
  channels =  new GXVDecayChannelVector;
}

GXDecayTable::~GXDecayTable()
{
  // remove and delete all contents  
  GXVDecayChannelVector::iterator iCh;
  for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
    delete (*iCh);
  }
  channels->clear();
  delete  channels;
  channels = nullptr;
  parent = nullptr;
}    

void GXDecayTable::Insert( GXVDecayChannel * aChannel){
  if (parent == nullptr) { 
    parent = (GXParticleDefinition*)(aChannel->GetParent()); 
  }
  if (parent != aChannel->GetParent()) {
#ifdef G4VERBOSE
    G4cout << " GXDecayTable::Insert :: bad GXVDecayChannel (mismatch parent) "
           << "       " << parent->GetParticleName()
           << " input:" << aChannel->GetParent()->GetParticleName() << G4endl;
#endif
  } else {
    G4double br = aChannel->GetBR();
    GXVDecayChannelVector::iterator iCh;
    for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
      if (br > (*iCh)->GetBR()) {
	channels->insert(iCh,aChannel);
	return;
      }
    }
    channels->push_back(aChannel);
  }
}

GXVDecayChannel *GXDecayTable::SelectADecayChannel(G4double parentMass)
{
  // check if contents exist
  if (channels->size()<1) return 0;

  if(parentMass<0.) parentMass=parent->GetPDGMass(); 

  GXVDecayChannelVector::iterator iCh;
  G4double sumBR = 0.;
  for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
    if ( !((*iCh)->IsOKWithParentMass(parentMass)) ) continue;
    sumBR += (*iCh)->GetBR();
  }
  if (sumBR <=0.0) {
#ifdef G4VERBOSE
    G4cout << " GXDecayTable::SelectADecayChannel :: no possible DecayChannel"
           << "       " << parent->GetParticleName() << G4endl;
#endif
    return nullptr;
  }

  const size_t MAX_LOOP = 10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    G4double sum = 0.0;
    G4double br= sumBR*G4UniformRand();
    // select decay channel
    for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
      sum += (*iCh)->GetBR();
      if ( !((*iCh)->IsOKWithParentMass(parentMass)) ) continue;
      if (br < sum) return (*iCh);
    }
  }
  return nullptr;
}

void GXDecayTable::DumpInfo() const
{
  G4cout << "GXDecayTable:  " << parent->GetParticleName() << G4endl;
  G4int index =0;
  GXVDecayChannelVector::iterator iCh;
  for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
    G4cout << index << ": ";
    (*iCh)->DumpInfo();
    index += 1;
  }
  G4cout << G4endl;
}











