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
// $Id: GXKaonZeroShort.cc 79342 2014-02-24 11:42:42Z gcosmo $
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

#include "GXKaonZeroShort.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"

#include "GXPhaseSpaceDecayChannel.hh"
#include "GXDecayTable.hh"

// ######################################################################
// ###                      KAONZEROSHORT                             ###
// ######################################################################

GXKaonZeroShort* GXKaonZeroShort::theInstance = 0;

GXKaonZeroShort* GXKaonZeroShort::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "kaon0S";
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
                 name,    0.497614*GeV,  7.3508e-12*MeV,       0.0,
                    0,              -1,             0,
                    1,               0,             0,
              "meson",               0,             0,         310,
                false,      0.08954*ns,          NULL,
                false,          "kaon",           310);

 //create Decay Table
  GXDecayTable* table = new GXDecayTable();

  // create decay channels
  GXVDecayChannel** mode = new GXVDecayChannel*[2];
  // kaon0s -> pi+ + pi-
  mode[0] = new GXPhaseSpaceDecayChannel("kaon0S",0.6920,2,"pi+","pi-");
  // kaon0s -> pi0 + pi0
  mode[1] = new GXPhaseSpaceDecayChannel("kaon0S",0.3069,2,"pi0","pi0");

  for (G4int index=0; index <2; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<GXKaonZeroShort*>(anInstance);
  return theInstance;
}

GXKaonZeroShort*  GXKaonZeroShort::KaonZeroShortDefinition()
{
  return Definition();
}

GXKaonZeroShort*  GXKaonZeroShort::KaonZeroShort()
{
  return Definition();
}

