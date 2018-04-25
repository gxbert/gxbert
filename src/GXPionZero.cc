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
// $Id: GXPionZero.cc 79357 2014-02-25 10:06:54Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "GXPionZero.hh"
#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"


#include "GXPhaseSpaceDecayChannel.hh"
//TK PionZero will not decay in BERT
//#include "G4DalitzDecayChannel.hh"
#include "GXDecayTable.hh"

// ######################################################################
// ###                          PIONZERO                              ###
// ######################################################################

GXPionZero* GXPionZero::theInstance = 0;

GXPionZero* GXPionZero::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "pi0";
  // search in particle table]
  GXParticleTable* pTable = GXParticleTable::GetParticleTable();
  GXParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding

   anInstance = new GXParticleDefinition(
                 name,   0.1349766*GeV,  7.73e-06*MeV,         0.0,
                    0,              -1,            +1,
                    2,               0,            -1,
              "meson",               0,             0,         111,
                false,      8.52e-8*ns,          NULL,
                false,            "pi",          111);

   // Life time is given from width
   anInstance->SetPDGLifeTime( hbar_Planck/(anInstance->GetPDGWidth()) );
     
//TK PionZero will not decay in BERT
/*
  //create Decay Table
  GXDecayTable* table = new GXDecayTable();

  // create a decay channel
  GXVDecayChannel* mode;
  // pi0 -> gamma + gamma
  mode = new GXPhaseSpaceDecayChannel("pi0",0.988,2,"gamma","gamma");
  table->Insert(mode);
  // pi0 -> gamma + e+ + e-
  mode = new G4DalitzDecayChannel("pi0",0.012,"e-","e+");
  table->Insert(mode);

   anInstance->SetDecayTable(table);
*/
  }
  theInstance = reinterpret_cast<GXPionZero*>(anInstance);
  return theInstance;
}

GXPionZero*  GXPionZero::PionZeroDefinition()
{
  return Definition();
}

GXPionZero*  GXPionZero::PionZero()
{
  return Definition();
}

