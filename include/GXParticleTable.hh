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
// $Id: GXParticleTable.hh 106143 2017-09-14 06:34:42Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	27 June 1996, H.Kurashige
// ------------------------------------------------------------
//      added fParticleMessenger         14 Nov., 97 H.Kurashige
//      added Create/DeleteMessenger     06 Jul., 98 H.Kurashige
//      modified FindIon                 02 Aug., 98 H.Kurashige
//      added dictionary for encoding    24 Sep., 98 H.Kurashige
//      added RemoveAllParticles()        8 Nov., 98 H.Kurashige
//         --------------------------------
//      fixed  some improper codings     08 Apr., 99 H.Kurashige
//      modified FindIon/GetIon methods  17 AUg., 99 H.Kurashige
//      implement new version for using STL map instaed of RW PtrHashedDictionary
//                                       28 ct., 99  H.Kurashige
//      modified implementation of Remove 21 Mar.,08  H.Kurashige
//      remove G4ShortLivedTable         25 July, 13 H.Kurashige
//      added support for MuonicAtom     September, 17 K.L.Genser
//

#ifndef GXParticleTable_h
#define GXParticleTable_h 1

#include <map>

#include "G4ios.hh"
#include "globals.hh"
//#include "G4Threading.hh"
#include "GXParticleDefinition.hh"
#include "GXParticleTableIterator.hh"

//class G4UImessenger;
//class G4ParticleMessenger;
class GXIonTable;

class GXParticleTable
{
 // Class Description
 //   GXParticleTable is the table of pointer to GXParticleDefinition
 //   GXParticleTable is a "singleton" (only one and staic object)
 //   In GXParticleTable, each GXParticleDefinition pointer is stored
 //   with its name as a key to itself. So, each  GXParticleDefinition 
 //   object must have unique name for itself. 
 //   

 public:

   typedef GXParticleTableIterator<G4String, GXParticleDefinition*>::Map G4PTblDictionary;
   typedef GXParticleTableIterator<G4String, GXParticleDefinition*> G4PTblDicIterator;
   typedef GXParticleTableIterator<G4int, GXParticleDefinition*>::Map G4PTblEncodingDictionary;
   typedef GXParticleTableIterator<G4int, GXParticleDefinition*> G4PTblEncodingDicIterator;

 protected:
   // default constructor
   GXParticleTable();
   // Copy constructor and assignment operator
   GXParticleTable(const  GXParticleTable &right);
   GXParticleTable & operator=(const GXParticleTable &);

 public:

   void SlaveGXParticleTable();
   void WorkerGXParticleTable();
   // This method is similar to the constructor. It is used by each worker
   // thread to achieve the partial effect as that of the master thread.

   virtual ~GXParticleTable();
   void DestroyWorkerGXParticleTable();
   // This method is similar to the destructor. It is used by each worker
   // thread to achieve the partial effect as that of the master thread.
  
 public: // With Description
   static GXParticleTable* GetParticleTable();
   // return the pointer to GXParticleTable object
   //   GXParticleTable is a "singleton" and can get its pointer by this function
   //   At the first time of calling this function, the GXParticleTable object
   //   is instantiated 

   G4bool   contains(const GXParticleDefinition *particle) const;
   G4bool   contains(const G4String &particle_name) const;
   // returns TRUE if the ParticleTable contains 

   G4int    entries() const;
   G4int    size() const;
   // returns the number of Particles in the ParticleTable
   
   GXParticleDefinition* GetParticle(G4int index) const;
   // returns a pointer to i-th particles in the ParticleTable
   //    0<= index < entries()

   const G4String& GetParticleName(G4int index) const;
   // returns name of i-th particles in the ParticleTable

   GXParticleDefinition* FindParticle(G4int  PDGEncoding );
   GXParticleDefinition* FindParticle(const G4String &particle_name);
   GXParticleDefinition* FindParticle(const GXParticleDefinition *particle);
   // returns a pointer to the particle (0 if not contained)

   GXParticleDefinition* FindAntiParticle(G4int  PDGEncoding );
   GXParticleDefinition* FindAntiParticle(const G4String &particle_name);
   GXParticleDefinition* FindAntiParticle(const GXParticleDefinition *particle);
   // returns a pointer to its anti-particle (0 if not contained)

   G4PTblDicIterator* GetIterator() const;
   // return the pointer of Iterator (RW compatible)
   
   void DumpTable(const G4String &particle_name = "ALL");
   // dump information of particles specified by name 
 
 public: //With Description

   GXIonTable* GetIonTable() const; 
   // return the pointer to G4IonTable object


 public: // With Description
   GXParticleDefinition* Insert(GXParticleDefinition *particle);
   // insert the particle into ParticleTable 
   // return value is same as particle if successfully inserted
   //              or pointer to another GXParticleDefinition object 
   //                   which has same name of particle
   //              or 0 if fail to insert by another reason

   GXParticleDefinition* Remove(GXParticleDefinition *particle);
   // Remove the particle from the table (not delete)

   void RemoveAllParticles();
   // remove all particles from GXParticleTable 

   void DeleteAllParticles();
   // remove and delete all particles from GXParticleTable  

 public:
//   G4UImessenger*       CreateMessenger();
//   void                 DeleteMessenger();
  // create/delete messenger for the particle table 
  // these methods are supposed to  be invoked by G4RunManager only
 
  protected:

   const G4PTblDictionary*    GetDictionary() const;

   const G4String& GetKey(const GXParticleDefinition *particle) const;
   // return key value of the particle (i.e. particle name)

   const G4PTblEncodingDictionary* GetEncodingDictionary() const; 
   // return the pointer to EncodingDictionary

 private:
   G4int verboseLevel;
   // controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

 public:
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;

//   static G4ThreadLocal G4ParticleMessenger* fParticleMessenger;
   static G4ThreadLocal G4PTblDictionary*  fDictionary;
   static G4ThreadLocal G4PTblDicIterator* fIterator;
   static G4ThreadLocal G4PTblEncodingDictionary* fEncodingDictionary;
   // These fields should be thread local or thread private. For a singleton
   // class, we can change any member field as static without any problem
   // because there is only one instance. Then we are allowed to add 
   // "G4ThreadLocal".
 
   //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
   //Phase I changes this member to be thread local 
   //,while each thread holds its own copy of particles.
   //Phase II changes this member back in order to share particles.
   static GXParticleTable*  fgParticleTable;

   static GXIonTable*            fIonTable;
   // This field should be thread private. However, we have to keep one copy
   // of the ion table pointer. So we change all important fields of G4IonTable
   // to the thread local variable.

   // These shadow pointers are used by each worker thread to copy the content
   // from the master thread. 

//   static G4ParticleMessenger* fParticleMessengerShadow;
   static G4PTblDictionary*  fDictionaryShadow;
   static G4PTblDicIterator* fIteratorShadow;
   static G4PTblEncodingDictionary* fEncodingDictionaryShadow;

 private:
   const G4String        noName;

   G4bool                readyToUse;
   GXParticleDefinition* genericIon;
   GXParticleDefinition* genericMuonicAtom;
 
 public:
   void SetReadiness(G4bool val=true);
   G4bool GetReadiness() const;
   GXParticleDefinition* GetGenericIon() const;
   void SetGenericIon(GXParticleDefinition*);
   GXParticleDefinition* GetGenericMuonicAtom() const;
   void SetGenericMuonicAtom(GXParticleDefinition*);
 private:
   void CheckReadiness() const;


#ifdef G4MULTITHREADED
public:
     //Andrea Dotti January 16. Shared instance of a mutex
     static G4Mutex particleTableMutex;
     static G4int lockCount;
#endif
};
#include "GXParticleTable.icc"

#endif
