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
// $Id: GXVShortLivedParticle.cc 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      History: first implementation, based on object model of
//      28 June 1998 H.Kurashige
// --------------------------------------------------------------

#include "GXVShortLivedParticle.hh"

GXVShortLivedParticle::GXVShortLivedParticle(const G4String&  aName,  
               G4double         mass,     
               G4double         width,
               G4double         charge,   
               G4int            iSpin,
               G4int            iParity,
               G4int            iConjugation,
               G4int            iIsospin,   
               G4int            iIsospinZ, 
               G4int            gParity,
               const G4String&  pType,
               G4int            lepton,
               G4int            baryon,
               G4int            encoding,
               G4bool           stable,
               G4double         lifetime,
	       GXDecayTable     *decaytable)
  :GXParticleDefinition( aName,mass,width,charge,iSpin,iParity,
           iConjugation,iIsospin,iIsospinZ,gParity,pType,
           lepton,baryon,encoding,stable,lifetime,decaytable, true)
{
   
}

GXVShortLivedParticle::~GXVShortLivedParticle() {}

G4int GXVShortLivedParticle::operator==(const GXVShortLivedParticle &right) const
{
  return (this->GetParticleName() == right.GetParticleName());
}

G4int GXVShortLivedParticle::operator!=(const GXVShortLivedParticle &right) const
{
  return (this->GetParticleName() != right.GetParticleName());
}