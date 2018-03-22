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
// $Id: GXDecayProducts.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      12 Dec 1997 H.Kurashige
//
//      Use std::vector 4  Apr. 2012  H.Kurashige
// ------------------------------------------------------------

#ifndef GXDecayProducts_h
#define GXDecayProducts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "GXDynamicParticle.hh"
#include  <vector>

class GXDecayProducts
{
  public: // With Description

    // constructors
    GXDecayProducts();
    GXDecayProducts(const GXDynamicParticle &aParticle);

  public: 
    // copy constructor and assignment operator 
    //   Deep    copy: for GXDynamicParticle pointer
    GXDecayProducts(const GXDecayProducts &right);
    GXDecayProducts & operator=(const GXDecayProducts &right);

    //destructor
    ~GXDecayProducts();

    // (un)equal operator
    G4int operator==(const GXDecayProducts &right) const;
    G4int operator!=(const GXDecayProducts &right) const;

  public: // With Description
   //  set-get methods for the parent particle   
   //    theParentPaticle is used to get information of parent particle 
   //    when decay products are filled 
   //    new GXDynamicParticle object is created in set methods  
    const GXDynamicParticle* GetParentParticle() const {return theParentParticle;};
    void SetParentParticle(const GXDynamicParticle &aParticle);

   //  boost all products
    void Boost(G4double totalEnergy, const G4ThreeVector &momentumDirection);
    void Boost(G4double betax, G4double betay, G4double betaz);
 
  //   push-pop  methods for decay products pointer
    GXDynamicParticle* PopProducts();
    G4int PushProducts(GXDynamicParticle *aParticle);

    GXDynamicParticle* operator[](G4int anIndex) const;

    G4int entries() const {return numberOfProducts;};

  // check energy/momentum of products 
    G4bool IsChecked() const; 
   
  // 
    void DumpInfo() const;

    typedef std::vector<GXDynamicParticle*>  G4DecayProductVector;
  protected:
    enum {MaxNumberOfProducts = 64};

  private: 
    G4int                               numberOfProducts;
    GXDynamicParticle*                  theParentParticle;
    G4DecayProductVector*               theProductVector;

};

// ------------------------
// Inlined operators
// ------------------------

inline 
 G4int GXDecayProducts::operator==(const GXDecayProducts &right) const
{
  return (this == (GXDecayProducts *) &right);
}

inline 
 G4int GXDecayProducts::operator!=(const GXDecayProducts &right) const
{
  return (this != (GXDecayProducts *) &right);
}

#endif

