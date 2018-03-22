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
// $Id: GXIonTable.hh 99621 2016-09-29 10:14:44Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, 
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      added clear()			20 Mar., 08 H.Kurashige
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      add GetNucleusMass              15 Mar. 99  H.Kurashige
//          -----
//      Modified GetIon methods                  17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Add GetNucleusEncoding according PDG 2006 9 Oct. 2006 H.Kurashige
//      Use STL map                              30 Jul. 2009 H.Kurashige
//      Add GetIsomerMass                       25 July 2013  H.Kurashige
//
#ifndef GXIonTable_h
#define GXIonTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "GXParticleDefinition.hh"
#include "GXParticleTable.hh"
#include "GXIons.hh"

#include <cmath>
#include <vector>
#include <map>

class GXParticleTable;
class G4VIsotopeTable; 
//class G4IsotopeProperty;
//class G4NuclideTable; 

class GXIonTable
{
 // Class Description
 //   GXIonTable is the table of pointer to GXParticleDefinition
 //   In GXIonTable, each GXParticleDefinition pointer is stored
 //

 public:
   // Use STL map as list of ions
   typedef  std::multimap<G4int, const GXParticleDefinition*> G4IonList;
   typedef  std::multimap<G4int, const GXParticleDefinition*>::iterator G4IonListIterator;

 public:
   //static GXIonTable* GetIonTable()
   //{ return G4ParticleTable::GetParticleTable()->GetIonTable(); }

 public:
  // constructor
   GXIonTable();

   void SlaveGXIonTable();
   void WorkerGXIonTable();
   // Method is used by each worker thread to copy the content from the master
   // thread.

 protected:
   // hide copy construictor as protected 
   GXIonTable(const  GXIonTable &right);
   GXIonTable & operator = (const GXIonTable &) {return *this;}

 public:
  // destructor
   virtual ~GXIonTable();
   void DestroyWorkerGXIonTable();

 public: // With Description
   G4int GetNumberOfElements() const;
   // Get number of elements defined in the IonTable

   // Register Isotope table
   // GXIonTable asks properties of isotopes to this G4VIsotopeTable 
   // by using FindIsotope(G4IsotopeProperty* property) method.
   
   // ---------------------------  
   // FindIon/GetIon
   //   FindIon methods return pointer of ion if it exists       
   //   GetIon methods also return pointer of ion. In GetIon 
   //   methods the designated ion will be created if it does not exist.
   //
   // !! PDGCharge inGXParticleDefinition of ions is           !!
   // !! electric charge of nucleus (i.e. fully ionized ions)  !!
   // -----------------------------

   // All ground state ions will be created
   //   stabele ground states are defined in G4NuclearProperty 
 
   // All excited ions with long life time (>1.0*ns) will be created
   //  isomers are defined in G4VIsotopeTable
   
   // All nuclide with a life time longer than certain value will be created
   // prior to the event loop.

   // Find/Get "ground state" and "excited state" 
   GXParticleDefinition* GetIon(G4int Z, G4int A, G4double E, G4int J=0);
   GXParticleDefinition* GetIon(G4int Z, G4int A, G4double E, 
                                GXIons::GXFloatLevelBase flb, G4int J=0);

   //   Z: Atomic Number
   //   A: Atomic Mass (nn + np +nlambda)
   //   L: Number of Lmabda
   //   E: Excitaion energy
   //   lvl:  Isomer Level 0: ground state)
   //   flb:  Floating level base (enum defined in GXIons.hh)
   //   flbChar:  Floating level base denoted by a character
   //             (<null>,X,Y,Z,U,V,W,R,S,T,A,B,C,D,E)
   //   J: Total Angular momentum (in unit of 1/2) : not used

   GXParticleDefinition* GetIon(G4int encoding);
   // The ion can be get by using PDG encoding 
   // !! Only ground state can be obtained .i.e. Isomer = 0
   
   // Find/Get "excited state" 
   GXParticleDefinition* FindIon(G4int Z, G4int A, G4double E, G4int J=0);
   GXParticleDefinition* FindIon(G4int Z, G4int A, G4double E,
                                 GXIons::GXFloatLevelBase flb, G4int J=0);

   //   Z: Atomic Number
   //   A: Atomic Mass (nn + np +nlambda)
   //   L: Number of Lmabda
   //   E: Excitaion energy
   //   lvl:  Isomer Level 0: ground state)
   //   flb:  Floating level base (enum defined in GXIons.hh)
   //   flbChar:  Floating level base denoted by a character
   //             (<null>,X,Y,Z,U,V,W,R,S,T,A,B,C,D,E)
   //   J: Total Angular momentum (in unit of 1/2) : not used

   static G4bool        IsIon(const GXParticleDefinition*);
   // return true if the particle is ion

   static G4bool        IsAntiIon(const GXParticleDefinition*);
   // return true if the particle is anti_ion


   const G4String&  GetIonName(G4int Z, G4int A, G4int lvl=0) const;
   const G4String&  GetIonName(G4int Z, G4int A, G4double E,
          GXIons::GXFloatLevelBase flb=GXIons::GXFloatLevelBase::no_Float) const;

   // get ion name
  
   static G4int GetNucleusEncoding(G4int Z,        G4int A, 
				   G4double E=0.0, G4int lvl=0);
  //  get PDG code for Ions 
  // Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn and Z = np.
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  
   static G4int GetNucleusEncoding(G4int Z,   G4int A,  G4int L,        
				   G4double E=0.0, G4int lvl=0);
  //  get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn +nlambda and Z = np.
  // L = nlambda
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations

   static G4bool GetNucleusByEncoding(G4int encoding,
				     G4int &Z,      G4int &A, 
				     G4double &E,   G4int &lvl);
   static G4bool GetNucleusByEncoding(G4int encoding,
				      G4int &Z,      G4int &A,  G4int &L,    
	 			      G4double &E,   G4int &lvl);
   // Energy will not be given even for excited state!!  
 
   G4double   GetNucleusMass(G4int Z, G4int A, G4int L=0, G4int lvl=0) const;
   // These methods returns Nucleus (i.e. full ionized atom) mass 
   // ,where Z is Atomic Number (number of protons) and
   //  A is Atomic Number (number of nucleons and hyperons)
   //  L is number of lambda (A= nn + np + nlambda)
   //  lvl is isomer level
 
   G4double   GetLifeTime(const GXParticleDefinition*) const;
   G4double   GetLifeTime(G4int Z, G4int A, G4double E,
                GXIons::GXFloatLevelBase flb=GXIons::GXFloatLevelBase::no_Float) const;
   G4double   GetLifeTime(G4int Z, G4int A, G4double E, char flbChar) const;
   // Returns a life time of an ion. -1 for stable ion, and -1001 for ion
   // that is not listed in G4NuclideTable.
   
   G4int                 Entries() const;
   // Return number of ions in the table

   GXParticleDefinition* GetParticle(G4int index) const;
   // Return the pointer of index-th ion in the table
 
   G4bool                Contains(const GXParticleDefinition *particle) const;
   // Return 'true' if the ion exists

   void                  Insert(const GXParticleDefinition* particle);
   void                  Remove(const GXParticleDefinition* particle);
   // Insert/Remove an ion in the table

   void                  clear();
   // erase all contents in the list (not delete just remove)

   G4int                 size() const;
   //  Return number of ions in the table

   void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 


 protected:
   GXParticleDefinition* FindIonInMaster(G4int Z, G4int A, G4int lvl=0);
   GXParticleDefinition* FindIonInMaster(G4int Z, G4int A, G4int L, G4int lvl);
   GXParticleDefinition* FindIonInMaster(G4int Z, G4int A, G4double E,
                                         GXIons::GXFloatLevelBase flb, G4int J=0);
   GXParticleDefinition* FindIonInMaster(G4int Z, G4int A, G4int L,
				 G4double E, GXIons::GXFloatLevelBase flb, G4int J=0);

   GXParticleDefinition* CreateIon(G4int Z, G4int A, G4double E, GXIons::GXFloatLevelBase flb);

   void                  InsertWorker(const GXParticleDefinition* particle);

   // Obsolete
   // GXParticleDefinition* SlaveCreateIon(GXParticleDefinition* ion, G4int Z, G4int A, G4double E);
   // All threads share the particle table and particles including ions. This method
   // is invoked by any work thread for ions that have been created by other threads
   // to achieve the partial effect when ions are created by other threads.
   // GXParticleDefinition* SlaveCreateIon(GXParticleDefinition* ion, G4int Z, G4int A, G4int L, G4double E);
   // All threads share the particle table and particles including ions. This method
   // is invoked by any work thread for ions that have been created by other threads
   // to achieve the partial effect when ions are created by other threads.

   // Create Ion 
   
   //TK for GXBERT
   //G4IsotopeProperty* FindIsotope(G4int Z, G4int A, G4double E, GXIons::GXFloatLevelBase flb) const;
   //G4IsotopeProperty* FindIsotope(G4int Z, G4int A, G4int  lvl) const; 
   // Ask properties of isotopes to this G4VIsotopeTable 
   
   GXParticleDefinition* GetLightIon(G4int Z, G4int A) const;
   GXParticleDefinition* GetLightAntiIon(G4int Z, G4int A) const;
   
   G4bool                IsLightIon(const GXParticleDefinition*) const;
   G4bool                IsLightAntiIon(const GXParticleDefinition*) const;
   // return true if the particle is pre-defined ion
 
   void                  AddProcessManager(GXParticleDefinition*);
   // Add process manager to ions with name of 'ionName'

   G4int                GetVerboseLevel() const;
   // get Verbose Level defined in G4ParticleTable

 private:
   //G4NuclideTable* pNuclideTable;
   G4bool         isIsomerCreated;
   // Isomer table and flag of creation    
 
 public:
   static G4ThreadLocal G4IonList* fIonList; 
   static G4ThreadLocal std::vector<G4VIsotopeTable*> *fIsotopeTableList;
   static G4IonList* fIonListShadow; 
   static std::vector<G4VIsotopeTable*> *fIsotopeTableListShadow;
   // It is very important for multithreaded Geant4 to keep only one copy of the
   // particle table pointer and the ion table pointer. However, we try to let 
   // each worker thread hold its own copy of the particle dictionary and the 
   // ion list. This implementation is equivalent to make the ion table thread
   // private. The two shadow ponters are used by each worker thread to copy the
   // content from the master thread.
 
   enum { numberOfElements = 118};
   static const G4String       elementName[numberOfElements];
    
   //needed for MT
   void InitializeLightIons();

 private:
   G4int n_error;

#ifdef G4MULTITHREADED
 public:
   static G4Mutex ionTableMutex;
#endif
};

inline G4int  GXIonTable::GetNumberOfElements() const
{
  return numberOfElements;
}

#endif
