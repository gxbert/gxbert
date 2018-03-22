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
// $Id: GXIons.hh 103892 2017-05-03 08:11:00Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 27 June 1998
// ----------------------------------------------------------------
//      Add excitation energy         17 Aug. 1999 H.Kurashige
//      Add isomer level              30 Apr. H.Kurashige


#ifndef GXIons_h
#define GXIons_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "GXParticleDefinition.hh"

// ######################################################################
// ###                          Ions                                 ###
// ######################################################################

class GXIons : public GXParticleDefinition
{
 // Class Description
 //  This is the base class for all nuclei including pre-defined 
 //  light nuclei such as deuteron, alpha, and proton (Hydrogen) 
 //  All nuclei/ions created on the fly are objects of this class
 //  Atomic number and atomic mass are vaild only for particles derived
 //  from this class.  This class has Excitation Energy in addition to
 //  the normal particle properties.

 protected:
   GXIons(){};


 public: //With Description
   GXIons(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       GXDecayTable        *decaytable,  G4bool              shortlived,
       const G4String&     subType ="",
       G4int               anti_encoding =0,
       G4double            excitation = 0.0, 
       G4int               isomer = 0
   );

 public:
   virtual    			~GXIons();
   GXIons*    			IonsDefinition();
   GXIons*    			Ions();

 public:  //With Description
   // Get excitation energy of nucleus
   G4double GetExcitationEnergy() const ; 
  
  // Get Isomer level (=0 for ground state)
  G4int GetIsomerLevel() const; 
   
  // enumerator for floating level base
  enum class GXFloatLevelBase
       { no_Float=0,
         plus_X, plus_Y, plus_Z, plus_U, plus_V, plus_W,
         plus_R, plus_S, plus_T, plus_A, plus_B, plus_C, plus_D, plus_E
       };
  static GXIons::GXFloatLevelBase FloatLevelBase(char flbChar);
  static GXIons::GXFloatLevelBase FloatLevelBase(G4int flbIdx);
  static char FloatLevelBaseChar(GXIons::GXFloatLevelBase flb);

  // set/get methods for floating level base
  GXIons::GXFloatLevelBase GetFloatLevelBase() const;
  G4int GetFloatLevelBaseIndex() const;
  void SetFloatLevelBase(GXIons::GXFloatLevelBase flb);
  void SetFloatLevelBase(char flbChar);
  void SetFloatLevelBase(G4int flbIdx);

 private:
  G4double theExcitationEnergy; 
  G4int    theIsomerLevel;
  GXFloatLevelBase floatLevelBase;

};

#define noFloat_ GXIons::GXFloatLevelBase::no_Float
#define plusU_ GXIons::GXFloatLevelBase::plus_U 
#define plusV_ GXIons::GXFloatLevelBase::plus_V 
#define plusW_ GXIons::GXFloatLevelBase::plus_W 
#define plusX_ GXIons::GXFloatLevelBase::plus_X
#define plusY_ GXIons::GXFloatLevelBase::plus_Y 
#define plusZ_ GXIons::GXFloatLevelBase::plus_Z 
#define plusR_ GXIons::GXFloatLevelBase::plus_R 
#define plusS_ GXIons::GXFloatLevelBase::plus_S 
#define plusT_ GXIons::GXFloatLevelBase::plus_T 
#define plusA_ GXIons::GXFloatLevelBase::plus_A
#define plusB_ GXIons::GXFloatLevelBase::plus_B 
#define plusC_ GXIons::GXFloatLevelBase::plus_C 
#define plusD_ GXIons::GXFloatLevelBase::plus_D 
#define plusE_ GXIons::GXFloatLevelBase::plus_E 

inline
 GXIons* GXIons::IonsDefinition()
{
  return this;
}

inline
 GXIons* GXIons::Ions() 
{
  return this;
}

inline
 G4double GXIons::GetExcitationEnergy() const 
{
  return theExcitationEnergy;
}

inline
 G4int GXIons::GetIsomerLevel() const
{
  return theIsomerLevel;
}
    
inline
 GXIons::GXFloatLevelBase GXIons::GetFloatLevelBase() const
{
  return floatLevelBase;
}

inline
 G4int GXIons::GetFloatLevelBaseIndex() const
{
  return static_cast<G4int>(floatLevelBase);
}

inline
 void GXIons::SetFloatLevelBase(GXIons::GXFloatLevelBase flb)
{
  floatLevelBase = flb;
}

inline
 void GXIons::SetFloatLevelBase(char flbChar)
{
  floatLevelBase = FloatLevelBase(flbChar);
}

inline
 void GXIons::SetFloatLevelBase(G4int flbIdx)
{
  floatLevelBase = FloatLevelBase(flbIdx);
}

inline
 GXIons::GXFloatLevelBase GXIons::FloatLevelBase(char flbChar)
{
  GXIons::GXFloatLevelBase flb = noFloat_;
  switch(flbChar)
  {
   case 'x': case 'X':
    flb = plusX_;
    break;
   case 'y': case 'Y':
    flb = plusY_;
    break;
   case 'z': case 'Z':
    flb = plusZ_;
    break;
   case 'u': case 'U':
    flb = plusU_;
    break;
   case 'v': case 'V':
    flb = plusV_;
    break;
   case 'w': case 'W':
    flb = plusW_;
    break;
   case 'r': case 'R':
    flb = plusR_;
    break;
   case 's': case 'S':
    flb = plusS_;
    break;
   case 't': case 'T':
    flb = plusT_;
    break;
   case 'a': case 'A':
    flb = plusA_;
    break;
   case 'b': case 'B':
    flb = plusB_;
    break;
   case 'c': case 'C':
    flb = plusC_;
    break;
   case 'd': case 'D':
    flb = plusD_;
    break;
   case 'e': case 'E':
    flb = plusE_;
    break;
   case '\0': default:
    break;
  }
  return flb;
}

inline
 GXIons::GXFloatLevelBase GXIons::FloatLevelBase(G4int flbIdx)
{
  static GXIons::GXFloatLevelBase flb[] = 
  { noFloat_,
    plusX_, plusY_, plusZ_, plusU_, plusV_, plusW_, 
    plusR_, plusS_, plusT_, plusA_, plusB_, plusC_, plusD_, plusE_ };
  return flb[flbIdx];
}

inline
 char GXIons::FloatLevelBaseChar(GXIons::GXFloatLevelBase flb)
{
  static char flbChar[] = {'\0','X','Y','Z','U','V','W',
                                'R','S','T','A','B','C','D','E'};
  return flbChar[static_cast<G4int>(flb)];
}

#endif








