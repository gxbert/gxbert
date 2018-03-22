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
// $Id: GXSigmaPlus.cc 67971 2013-03-13 10:13:24Z gcosmo $
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

#include "GXSigmaPlus.hh"
#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"

#include "GXPhaseSpaceDecayChannel.hh"
#include "GXDecayTable.hh"

// ######################################################################
// ###                           SigmaPlus                            ###
// ######################################################################

GXSigmaPlus* GXSigmaPlus::theInstance = 0;

GXSigmaPlus* GXSigmaPlus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "sigma+";
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
                 name,    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
                false,      0.08018*ns,          NULL,
                false,       "sigma");
 
    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( 2.458 * mN);

    //create Decay Table
    GXDecayTable* table = new GXDecayTable();
    
    // create decay channels
    GXVDecayChannel** mode = new GXVDecayChannel*[2];
    // sigma+ -> proton + pi0
    mode[0] = new GXPhaseSpaceDecayChannel("sigma+",0.516,2,"proton","pi0");
    // sigma+ -> neutron + pi+
    mode[1] = new GXPhaseSpaceDecayChannel("sigma+",0.483,2,"neutron","pi+");
    
    for (G4int index=0; index <2; index++ ) table->Insert(mode[index]);
    delete [] mode;
    
    anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<GXSigmaPlus*>(anInstance);
  return theInstance;
}

GXSigmaPlus*  GXSigmaPlus::SigmaPlusDefinition()
{
  return Definition();
}

GXSigmaPlus*  GXSigmaPlus::SigmaPlus()
{
  return Definition();
}

