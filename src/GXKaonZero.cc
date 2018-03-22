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
// $Id: GXKaonZero.cc 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
//                              H.Kurashige   7 Jul 96
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "GXKaonZero.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"

#include "GXPhaseSpaceDecayChannel.hh"
#include "GXDecayTable.hh"

// ######################################################################
// ###                          KAONZERO                              ###
// ######################################################################

GXKaonZero* GXKaonZero::theInstance = 0;

GXKaonZero* GXKaonZero::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "kaon0";
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
                 name,    0.497614*GeV,       0.0*MeV,         0.0,
                    0,              -1,             0,
                    1,              -1,             0,
              "meson",               0,             0,         311,
                false,             0.0,          NULL,
                false,       "kaon");

 //create Decay Table
  GXDecayTable* table = new GXDecayTable();

 // create decay channels
  GXVDecayChannel** mode = new GXVDecayChannel*[2];
  // kaon0 -> Kaon0L
  mode[0] = new GXPhaseSpaceDecayChannel("kaon0",0.500,1,"kaon0L");
  // kaon0 -> Kaon0S
  mode[1] = new GXPhaseSpaceDecayChannel("kaon0",0.500,1,"kaon0S");

  for (G4int index=0; index <2; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<GXKaonZero*>(anInstance);
  return theInstance;
}

GXKaonZero*  GXKaonZero::KaonZeroDefinition()
{
  return Definition();
}

GXKaonZero*  GXKaonZero::KaonZero()
{
  return Definition();
}

