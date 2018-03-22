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
// $Id: GXDecayTable.hh 105720 2017-08-16 12:38:10Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      7 July 1996 H.Kurashige
//
// ----------------------------------------
//      implementation for STL          14 Feb. 2000 H.Kurashige
//      
// ------------------------------------------------------------

#ifndef GXDecayTable_h
#define GXDecayTable_h 1

#include "G4ios.hh"
#include <vector>
#include "globals.hh"
#include "GXParticleDefinition.hh"
#include "GXVDecayChannel.hh"

class GXDecayTable 
{
 // Class Description
 //   GXDecayTable is the table of pointer to GXVDecayChannel.
 //   Decay channels inside is sorted by using decay branching ratio
 //

 public:
  typedef std::vector<GXVDecayChannel*> GXVDecayChannelVector;

  //constructors
 public:
    GXDecayTable();
    ~GXDecayTable();

 private:
  // hide copy constructor and assignment operator by declaring "private"
  //  (Implementation does not make sense )
    GXDecayTable(const GXDecayTable &){};
    GXDecayTable & operator=(const GXDecayTable &){return *this;};

 public:
    // equality operators
    G4int operator==(const GXDecayTable &right) const {return (this == &right);};
    G4int operator!=(const GXDecayTable &right) const {return (this != &right);};

 public: // With Description
    void  Insert( GXVDecayChannel* aChannel);
    // Insert a decay channel at proper position 
    // (i.e. sorted by using branching ratio ) 

    G4int entries() const;
    // Returns number of decay channels inside

 public: // With Description
    GXVDecayChannel* SelectADecayChannel(G4double parentMass= -1.);
    // A decay channel is selected at random according to the branching ratio 
  
    GXVDecayChannel* GetDecayChannel(G4int index) const;
    GXVDecayChannel* operator[](G4int index);
    // Get index-th Decay channel
 
    void DumpInfo() const;

 private:
    GXParticleDefinition       *parent;
    GXVDecayChannelVector      *channels;
};

inline     
 G4int GXDecayTable::entries() const
{
  return channels->size();
}

inline     
 GXVDecayChannel* GXDecayTable::operator[](G4int index)
{
  return (*channels)[index];
}

 
inline     
 GXVDecayChannel* GXDecayTable::GetDecayChannel(G4int index) const
{
  GXVDecayChannel* selectedChannel = nullptr;
  if ( (index>=0) && (index<G4int(channels->size())) ){
    selectedChannel = (*channels)[index];
  }
  return selectedChannel;
}
 
 
#endif
