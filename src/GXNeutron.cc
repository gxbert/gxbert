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
// $Id: GXNeutron.cc 102905 2017-03-02 09:50:56Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
//                          H.Kurashige 7 July 1996
//      add neutron life time    Oct 17 2000 
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "GXNeutron.hh"
#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"

//TK neutron will not decay in BERT
//#include "GXNeutronBetaDecayChannel.hh"
#include "GXDecayTable.hh"

// ######################################################################
// ###                           NEUTRON                              ###
// ######################################################################
GXNeutron* GXNeutron::theInstance = 0;

GXNeutron* GXNeutron::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "neutron";
  // search in particle table]
  GXParticleTable* pTable = GXParticleTable::GetParticleTable();
  GXIons* anInstance =  reinterpret_cast<GXIons*>(pTable->FindParticle(name));
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
  // use constants in CLHEP
  // static const double  neutron_mass_c2 = 939.56563 * MeV;

    anInstance = new GXIons(
                 name, neutron_mass_c2, 7.478e-28*GeV,         0.0, 
		    1,              +1,             0,          
		    1,              -1,             0,             
	     "baryon",               0,            +1,        2112,
		false,    880.2*second,          NULL,
		false,       "nucleon",         -2112,
		 0.0,                0 
              );
    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( -1.9130427 * mN);
    //create Decay Table 
//TK neutron will not decay in BERT
/*
    GXDecayTable* table = new GXDecayTable();
    // create a decay channel
    GXVDecayChannel* mode = new GXNeutronBetaDecayChannel("neutron",1.00);
    table->Insert(mode);
    anInstance->SetDecayTable(table);
*/
    
  }
  theInstance = reinterpret_cast<GXNeutron*>(anInstance);
  return theInstance;
}

GXNeutron*  GXNeutron::NeutronDefinition()
{
  return Definition();
}

GXNeutron*  GXNeutron::Neutron()
{
  return Definition();
}

