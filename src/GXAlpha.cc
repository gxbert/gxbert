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
// $Id: GXAlpha.cc 69557 2013-05-08 12:01:40Z gcosmo $
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

#include "GXAlpha.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"

// ######################################################################
// ###                           ALPHA                                ###
// ######################################################################

GXAlpha* GXAlpha::theInstance = 0;

GXAlpha* GXAlpha::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "alpha";
  // search in particle table]
  GXParticleTable* pTable = GXParticleTable::GetParticleTable();
  GXIons* anInstance = reinterpret_cast<GXIons*>(pTable->FindParticle(name));
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
  //             excitation
   anInstance = new GXIons(
                 name,    3.727379*GeV,       0.0*MeV,  +2.0*eplus,
                    0,              +1,             0,
                    0,               0,             0,
            "nucleus",               0,            +4,  1000020040,
                 true,            -1.0,          NULL,
		false,        "static",   -1000020040, 
                  0.0,               0
               );

  }

  theInstance = reinterpret_cast<GXAlpha*>(anInstance);
  return theInstance;
}

GXAlpha*  GXAlpha::AlphaDefinition()
{
  return Definition();
}

GXAlpha*  GXAlpha::Alpha()
{
  return Definition();
}


