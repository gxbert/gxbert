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
// $Id: GXElectronOccupancy.hh 79357 2014-02-25 10:06:54Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 17 Aug 1999
// ----------------------------------------------------------------
// Class Description
//     This class has information of occupation of electrons 
//     in atomic orbits
// -  
//     GetOccupancy(N) gives the number of electron
//     in N-th orbit
//       For example : Carbon atom should be 
//          GetOccupancy(0)      --> 2
//          GetOccupancy(1)      --> 4
//          GetOccupancy(2..7)   --> 0
// -
//     GetTotalOccupancy() gives the total number of electrons
//
// --------------------------------------------------------------- 



#ifndef GXElectronOccupancy_h
#define GXElectronOccupancy_h 1

#include "globals.hh"
#include "GXAllocator.hh"
#include "G4ios.hh"

//#include "pwdefs.hh"
#define G4PART_DLL

class GXElectronOccupancy 
{
 public:
   enum { MaxSizeOfOrbit = 20};

 public: // With Description
   GXElectronOccupancy( G4int sizeOrbit = MaxSizeOfOrbit   );
   GXElectronOccupancy( const GXElectronOccupancy& right );

 public:
   virtual    	       ~GXElectronOccupancy();

  //  new/delete operators are oberloded to use GXAllocator
     inline void *operator new(size_t);
     inline void operator delete(void *aGXElectronOccupancy);

 
  //- operators
     GXElectronOccupancy & operator=(const GXElectronOccupancy &right);
     G4int operator==(const GXElectronOccupancy &right) const;
     G4int operator!=(const GXElectronOccupancy &right) const;
   
 public: // With Description
   // The following methods returns
   //     0:  if the orbit(atom) is vacant 
   //    >0:  number of electrons in orbit
   G4int  GetTotalOccupancy() const;
   G4int  GetOccupancy(G4int orbit) const;
 
   //
   G4int  AddElectron(G4int orbit, G4int number = 1);
   G4int  RemoveElectron(G4int orbit, G4int number = 1);
   
   G4int  GetSizeOfOrbit() const;
   void   DumpInfo() const;

 private:
   G4int  theSizeOfOrbit;
   G4int  theTotalOccupancy;
   G4int* theOccupancies;

};

extern G4PART_DLL G4ThreadLocal GXAllocator<GXElectronOccupancy> *aGXElectronOccupancyAllocator;

// ------------------------
// Inlined operators
// ------------------------

inline void * GXElectronOccupancy::operator new(size_t)
{
  if (!aGXElectronOccupancyAllocator)
  {
    aGXElectronOccupancyAllocator = new GXAllocator<GXElectronOccupancy>;
  }
  return (void *) aGXElectronOccupancyAllocator->MallocSingle();
}

inline void GXElectronOccupancy::operator delete(void * aGXElectronOccupancy)
{
  aGXElectronOccupancyAllocator->FreeSingle((GXElectronOccupancy *) aGXElectronOccupancy);
}

inline
 G4int  GXElectronOccupancy::GetSizeOfOrbit() const
{
  return  theSizeOfOrbit;
}

inline
 G4int GXElectronOccupancy::GetTotalOccupancy() const
{
  return  theTotalOccupancy;
}

inline
 G4int  GXElectronOccupancy::GetOccupancy(G4int orbit) const
{
  G4int value = 0;
  if ((orbit >=0)&&(orbit<theSizeOfOrbit)){
    value = theOccupancies[orbit];
  }
  return value;  
}


#endif
