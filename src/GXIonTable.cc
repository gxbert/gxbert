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
// $Id: GXIonTable.cc 106143 2017-09-14 06:34:42Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	27 June 1998  H.Kurashige
// ---------------------------------------------------------------
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      use G4NucleiPropoerties to get nuceli Mass 17  Nov.,98 H.Kurashige
//      use G4GenericIon for process List
//      modify fomula of Ion mass       09 Dec., 98 H.Kurashige 
//          -----
//      Modified GetIon methods         17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Modified Element Name for Z>103  06 Apr. 01 H.Kurashige
//      Remove test of cuts in SetCuts   16 Jan  03 V.Ivanchenko
//      Added initial support for Muonic Atoms   1 Jul 16  K.Lynch
//      Extended support for Muonic Atoms        September 17  K.L.Genser

#include <iostream>               
#include <iomanip>               
#include <sstream>

#include "G4ios.hh"
//#include "G4Threading.hh"

#include "GXIonTable.hh"
#include "GXPhysicalConstants.hh"
#include "GXSystemOfUnits.hh"
#include "GXParticleTable.hh"
//#include "G4StateManager.hh"
#include "GXIons.hh"
#include "GXNucleiProperties.hh"
//#include "G4HyperNucleiProperties.hh"

//#include "G4IsotopeProperty.hh"
//#include "G4VIsotopeTable.hh"
//#include "G4NuclideTable.hh"

// It is very important for multithreaded Geant4 to keep only one copy of the
// particle table pointer and the ion table pointer. However, we try to let 
// each worker thread hold its own copy of the particle dictionary and the 
// ion list. This implementation is equivalent to make the ion table thread
// private. The two shadow ponters are used by each worker thread to copy the
// content from the master thread.
//
G4ThreadLocal GXIonTable::G4IonList* GXIonTable::fIonList = 0;
G4ThreadLocal std::vector<G4VIsotopeTable*> *GXIonTable::fIsotopeTableList = 0;
GXIonTable::G4IonList* GXIonTable::fIonListShadow = 0;
std::vector<G4VIsotopeTable*> *GXIonTable::fIsotopeTableListShadow = 0;

namespace lightions {
  static const GXParticleDefinition* p_proton=0;
  static const GXParticleDefinition* p_deuteron=0;
  static const GXParticleDefinition* p_triton=0;
  static const GXParticleDefinition* p_alpha=0;
  static const GXParticleDefinition* p_He3=0;
  void Init() {
    if ( p_proton ) return;
    p_proton   = GXParticleTable::GetParticleTable()-> FindParticle("proton"); // proton
    p_deuteron = GXParticleTable::GetParticleTable()-> FindParticle("deuteron"); // deuteron
    p_triton   = GXParticleTable::GetParticleTable()-> FindParticle("triton"); // tritoon
    p_alpha    = GXParticleTable::GetParticleTable()-> FindParticle("alpha"); // alpha
    p_He3      = GXParticleTable::GetParticleTable()-> FindParticle("He3"); // He3
  }
}

namespace antilightions {
    static const GXParticleDefinition* p_proton=0;
    static const GXParticleDefinition* p_deuteron=0;
    static const GXParticleDefinition* p_triton=0;
    static const GXParticleDefinition* p_alpha=0;
    static const GXParticleDefinition* p_He3=0;
    void Init() {
        if ( p_proton ) return;
        p_proton   = GXParticleTable::GetParticleTable()-> FindParticle("anti_proton"); // proton
        p_deuteron = GXParticleTable::GetParticleTable()-> FindParticle("anti_deuteron"); // deuteron
        p_triton   = GXParticleTable::GetParticleTable()-> FindParticle("anti_triton"); // tritoon
        p_alpha    = GXParticleTable::GetParticleTable()-> FindParticle("anti_alpha"); // alpha
        p_He3      = GXParticleTable::GetParticleTable()-> FindParticle("anti_He3"); // He3
    }
}

////////////////////
GXIonTable::GXIonTable()
  : //pNuclideTable(0),
    isIsomerCreated(false),
    n_error(0)
{
  fIonList = new G4IonList();

  // Set up the shadow pointer used by worker threads.
  //
  if (fIonListShadow == 0)
  {
    fIonListShadow = fIonList;
  }

  fIsotopeTableList = new std::vector<G4VIsotopeTable*>;

  // Set up the shadow pointer used by worker threads.
  //
  if (fIsotopeTableListShadow == 0)
  {
    fIsotopeTableListShadow = fIsotopeTableList;
  }    

  //PrepareNuclideTable();
  //RegisterIsotopeTable(pNuclideTable);
}

// This method is used by each worker thread to copy the content
// from the master thread.
//
void GXIonTable::SlaveGXIonTable()
{
G4Exception("GXIonTable::SlaveGXParticleTable()","G4MT0000",FatalException,"Obsolete");
}

void GXIonTable::WorkerGXIonTable()
{
  if( fIonList == 0 )
  { fIonList = new G4IonList(); }
  else
  { fIonList->clear(); }

  G4IonListIterator it;
  for (it = fIonListShadow->begin() ; it != fIonListShadow->end(); it++ ) {
///////////////////////////    GXParticleDefinition* ion = const_cast<GXParticleDefinition*>(it->second);
///////////////////////////    if (ion->IsGeneralIon()) AddProcessManager(ion); 
    fIonList->insert(*it);
  }

  // Do not copy Isotoper Table to Worker thread
  if( fIsotopeTableList == 0 ) {
    fIsotopeTableList = new std::vector<G4VIsotopeTable*>; 
    for (size_t i = 0; i < fIsotopeTableListShadow->size(); i++){
      fIsotopeTableList->push_back((*fIsotopeTableListShadow)[i]); 
    }
  }

  /////////fIsotopeTableList = new std::vector<G4VIsotopeTable*>;
  /////////RegisterIsotopeTable(pNuclideTable);
}

void GXIonTable::InitializeLightIons()
{
  lightions::Init();
  antilightions::Init();
}


////////////////////
GXIonTable::~GXIonTable()
{
  // delete IsotopeTable if exists
  if (fIsotopeTableList != 0)
  {
    for (size_t i = 0; i< fIsotopeTableList->size(); ++i)
    {
      G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
      //delete fIsotopeTable;
      //TK for GXBERT   
      //if( fIsotopeTable != G4NuclideTable::GetNuclideTable() ) delete fIsotopeTable;
      if( fIsotopeTable != NULL ) delete fIsotopeTable;
    }
    fIsotopeTableList->clear();
    delete fIsotopeTableList;
  }
  fIsotopeTableList =0;


  if (fIonList ==0) return;
  // remove all contents in the Ion List 
  //  No need to delete here because all particles are dynamic objects
  fIonList->clear();

  delete fIonList;
  fIonList =0;
}

////////////////////
void GXIonTable::DestroyWorkerGXIonTable()
{
  // delete IsotopeTable if exists
  if (fIsotopeTableList != 0)
  {
    for (size_t i = 0; i< fIsotopeTableList->size(); ++i)
    {
      G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
      //delete fIsotopeTable;
      //TK for GXBERT   
      //if( fIsotopeTable != G4NuclideTable::GetNuclideTable() ) delete fIsotopeTable;
      if( fIsotopeTable != NULL ) delete fIsotopeTable;
    }
    fIsotopeTableList->clear();
    delete fIsotopeTableList;
  }
  fIsotopeTableList =0;


  if (fIonList ==0) return;
  // remove all contents in the Ion List 
  //  No need to delete here because all particles are dynamic objects
  fIonList->clear();

  delete fIonList;
  fIonList =0;
}


////////////////////
// -- CreateIon method ------
////////////////////
GXParticleDefinition* GXIonTable::CreateIon(G4int Z, G4int A, G4double E, 
                                            GXIons::GXFloatLevelBase flb)
{
  GXParticleDefinition* ion=0;

  // check whether GenericIon has processes
  GXParticleDefinition* genericIon = 
    GXParticleTable::GetParticleTable()->GetGenericIon();
//TK GXIons is internal BERT, following check is not necessary
/*
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (genericIon->GetParticleDefinitionID() < 0) || (pman==0)){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "GXIonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A 
             << "  because GenericIon is not ready !!" <<   G4endl;
    }
#endif
    G4Exception( "GXIonTable::CreateIon()","PART105",
		 JustWarning, 
		 "Can not create ions because GenericIon is not ready");
    return 0;
  }
*/
  
  G4double life = 0.0;
  GXDecayTable* decayTable =0;
  G4bool stable = true;
  G4double mu = 0.0;
  G4double Eex = 0.0;
  G4int    lvl =0;
  G4int    J=0;

//TK for GXBERT
  Eex = E;
  if (Eex>0.0) lvl=9;
/*
  const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, E, flb);
  if (fProperty !=0 ){
    Eex  = fProperty->GetEnergy();
    lvl  = fProperty->GetIsomerLevel();
    J    = fProperty->GetiSpin();
    life = fProperty->GetLifeTime();
    mu   = fProperty->GetMagneticMoment();    
    //TK do not consider decay of ions 
    //decayTable = fProperty->GetDecayTable();
    decayTable = NULL;
    stable = (life <= 0.) || (decayTable ==0);
    lvl = fProperty->GetIsomerLevel();
    if (lvl <0) lvl=9;
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4ExceptionDescription ed;
      ed << "GXIonTable::CreateIon() : G4IsotopeProperty object was not found for"
         << " Z = " << Z << " A = " << A << " E = " << E/keV << " (keV)";
      if(flb!=GXIons::GXFloatLevelBase::no_Float)
      { ed << " FloatingLevel +" << GXIons::FloatLevelBaseChar(flb); }
      ed << ".\n"
         << " Physics quantities such as life are not set for this ion.";
      G4Exception( "GXIonTable::CreateIon()","PART70105", JustWarning, ed);
    }
#endif
    // excitation energy
    Eex = E;
    // lvl is assigned to 9 temporally    
    if (Eex>0.0) lvl=9;
  }
*/

  //Eex = G4NuclideTable::Round(Eex); 
  if (Eex==0.0) lvl=0;
  // ion name
  G4String name =""; 
  /////////////if (lvl<9) name = GetIonName(Z, A, lvl);
  if (lvl==0 && flb==GXIons::GXFloatLevelBase::no_Float) name = GetIonName(Z, A, lvl);
  else       name = GetIonName(Z, A, Eex, flb);

  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,E,lvl);

//G4cout<<"GXIonTable::CreateIon "<<"Z:"<<Z<<" A:"<<A<<" E:"<<E<<" Eex:"<<Eex<<" lvl:"<<lvl<<" name:"<<name<<" code:"<<encoding<<G4endl;
  // PDG mass
  G4double mass =  GetNucleusMass(Z, A)+ Eex;
 
  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;
 
  // create an ion
  //   spin, parity, isospin values are fixed

  // Request lock for particle table accesses. Some changes are inside 
  // this critical region.
  //

  ion = new GXIons(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		 "generic",               0,
		       Eex,             lvl         );

  // Release lock for particle table accesses.
  //

  ion->SetPDGMagneticMoment(mu);
  static_cast<GXIons*>(ion)->SetFloatLevelBase(flb);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "GXIonTable::CreateIon() : create ion of " << name
	   << "  " << Z << ", " << A
	   << " encoding=" << encoding;
    if (E>0.0) {
      G4cout << " IsomerLVL=" << lvl
	     << " excited energy=" << Eex/keV << "[keV]";
    }
    G4cout << G4endl;
  } 
#endif
  
  // Add process manager to the ion
  //TK GXIon is an interanl BERT object, thus it does not need process manager
  //AddProcessManager(ion);
 
  return ion;
}


////////////////////

////////////////////
GXParticleDefinition* GXIonTable::GetIon(G4int Z, G4int A, G4double E, G4int J)
{ return GetIon(Z,A,E,GXIons::GXFloatLevelBase::no_Float,J); }

////////////////////
GXParticleDefinition* GXIonTable::GetIon(G4int Z, G4int A, G4double E,
                          GXIons::GXFloatLevelBase flb, G4int J)
{
  if ( (A<1) || (Z<=0) || (E<0.0) || (A>999) || (J<0) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "GXIonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  GXParticleDefinition* ion = FindIon(Z,A,E,flb,J);

  // create ion
  if (ion == 0) ion = CreateIon(Z,A,E,flb);

  return ion;  
}

////////////////////
GXParticleDefinition* GXIonTable::GetIon(G4int encoding)
{
  G4int Z, A, LL, IsoLvl;
  G4double E;
  if (!GetNucleusByEncoding(encoding,Z,A,LL,E,IsoLvl) ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "GXIonTable::GetIon() : illegal encoding" 
             << " CODE:" << encoding << G4endl;
    }
#endif
    G4Exception( "GXIonTable::GetIon()","PART106",
		 JustWarning, "illegal encoding for an ion");
    return 0;
  }
  //
  return GetIon( Z, A, LL, IsoLvl);
}

/////////////////////
// -- FindIon methods  ------
/////////////////////
GXParticleDefinition* GXIonTable::FindIon(G4int Z, G4int A, G4double E, G4int J)
{ return FindIon(Z,A,E,GXIons::GXFloatLevelBase::no_Float,J); }

////////////////////
GXParticleDefinition* GXIonTable::FindIon(G4int Z, G4int A, G4double E,
                           GXIons::GXFloatLevelBase flb, G4int J)
{
  if ( (A<1) || (Z<=0) || (J<0) || (E<0.0) || (A>999) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "GXIonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    G4Exception( "GXIonTable::FindIon()","PART107",
		 JustWarning, "illegal atomic number/mass");
    return 0;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const GXParticleDefinition* ion=0;
  G4bool isFound = false;

  // check if light ion
  ion = GetLightIon(Z,A);
  if (ion!=0 && E==0.0) { 
    // light ion 
    isFound = true;
  } else {    
    // -- loop over all particles in Ion table
    G4int encoding=GetNucleusEncoding(Z, A);
    G4IonList::iterator i = fIonList->find(encoding);
    for( ;i != fIonList->end() ; i++) {
      ion = i->second;
      if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
      // excitation level
      G4double anExcitaionEnergy = ((const GXIons*)(ion))->GetExcitationEnergy();
      //TK for GXBERT
      //if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
      if (std::fabs(E - anExcitaionEnergy) < 1.0*CLHEP::eV ) {
        if(((const GXIons*)(ion))->GetFloatLevelBase()==flb)
        {
	  isFound = true;
	  break;
        }
      }
    }
  }

  if ( isFound ){ 
    return const_cast<GXParticleDefinition*>(ion);
  } else {
    return 0;
  }
}

/////////////////
G4int GXIonTable::GetNucleusEncoding(G4int Z, G4int A, G4double E, G4int lvl)
{
  // PDG code for Ions
  // Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn and Z = np.
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  
  if ( Z==1 && A==1 && E==0.0 ) return 2212; // proton
  
  G4int encoding = 1000000000;
  encoding += Z * 10000;
  encoding += A *10;
  if (lvl>0&&lvl<10) encoding +=lvl;       //isomer level
  else if (E>0.0) encoding += 9;  //isomer level
  
  return encoding;
}

/////////////////
G4int GXIonTable::GetNucleusEncoding(G4int Z,  G4int A,    G4int LL,
				     G4double E,    G4int lvl)
{
  //  get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn +nlambda and Z = np.
  // LL = nlambda
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations

  G4int encoding = GetNucleusEncoding(Z, A, E, lvl);
  if (LL==0) return encoding; 
  encoding += LL*  10000000;
  if ( Z==1 && A==1 && E==0.0 ) encoding = 3122; // Lambda 

  return encoding;
}

///////////////
G4bool GXIonTable::GetNucleusByEncoding(G4int encoding,
			    G4int &Z,      G4int &A, 
			    G4double &E,   G4int &lvl)
{
  if (encoding <= 0) return false; // anti particle   

  if (encoding == 2212) {  // proton
    Z = 1; A = 1;
    E = 0.0; lvl =0; 
    return true;
  } 

  encoding -= 1000000000;
  Z = encoding/10000;
  encoding -= 10000*Z;
  A = encoding/10;
  lvl = encoding % 10;
  return true;
}

///////////////
G4bool GXIonTable::GetNucleusByEncoding(G4int encoding,
					G4int &Z,      G4int &A, 
					G4int &LL,   
					G4double &E,   G4int &lvl)
{
  if (encoding <= 0) return false; // anti particle   

 if (encoding == 3122) {  // Lambda
   Z = 1; A = 1; LL = 1;
    E = 0.0; lvl =0; 
    return true;
  } 

  if (encoding % 10 != 0) {
    //!!!not supported for excitation states !!!   
    return false;
  }
  if (encoding < 1000000000) {
    // anti particle   
    return false;
  }

  encoding -= 1000000000;
  LL = encoding/10000000;
  encoding -= 10000000*LL;
  Z = encoding/10000;
  encoding -= 10000*Z;
  A = encoding/10;
  lvl = encoding % 10;
  return true;
}

//TK for GXBERT: Exclude to use G4AutoDelete
//#include "G4AutoDelete.hh"
/////////////////
const G4String& GXIonTable::GetIonName(G4int Z, G4int A, G4double E,
                GXIons::GXFloatLevelBase flb) const 
{
  static G4ThreadLocal G4String *pname = 0;
                                       //TK for GXBERT: Exclude to use G4AutoDelete
  if (!pname)  { pname = new G4String("");/* G4AutoDelete::Register(pname);*/ }
  G4String &name = *pname;

  static G4ThreadLocal std::ostringstream* os = 0;
  if ( ! os ) {
    os = new std::ostringstream();
    //TK for GXBERT: Exclude to use G4AutoDelete
    //G4AutoDelete::Register(os); 
    os->setf(std::ios::fixed);
    os->precision(3);
  }

  name = GetIonName(Z, A);

  //excited energy
  if ( E>0  || flb!=GXIons::GXFloatLevelBase::no_Float){
    os->str("");
    std::ostringstream& oo = *os;
    // Excited nucleus
    oo<<'['<<E/keV;
    if(flb!=GXIons::GXFloatLevelBase::no_Float)
    { oo<<GXIons::FloatLevelBaseChar(flb); }
    oo<< ']';
    name += os->str();
  }

  return name;
}

/////////////////
const G4String& GXIonTable::GetIonName(G4int Z, G4int A, G4int lvl) const 
{
  static G4ThreadLocal G4String *pname = 0;
                                         //TK for GXBERT: Exclude to use G4AutoDelete
  if (!pname)  { pname = new G4String("");/* G4AutoDelete::Register(pname);*/ }
  G4String &name = *pname;

  static G4ThreadLocal std::ostringstream* os = 0;
  if ( ! os ) {
    os = new std::ostringstream();
    //TK for GXBERT: Exclude to use G4AutoDelete
    //G4AutoDelete::Register(os);
    os->setf(std::ios::fixed);
  }

  if ( (0< Z) && (Z <=numberOfElements) ) {
    name = elementName[Z-1];
  } else if (Z > numberOfElements) {
    os->str("");
    os->operator<<(Z);
    name = "E" + os->str() + "-";
  } else {
    name = "?";
    return name;
  }
  // Atomic Mass
  os->str("");
  os->operator<<(A);

  if ( lvl>0 ){
    std::ostringstream& oo = *os;
    // isomer level for Excited nucelus 
    oo<<'['<<lvl << ']';
  }
  name += os->str();

  return name;
}

/////////////////
G4bool GXIonTable::IsIon(const GXParticleDefinition* particle)
{
  // return true if the particle is ion

  static const G4String nucleus("nucleus");
  static const G4String proton("proton");

  // neutron is not ion
  if ((particle->GetAtomicMass()>0)   && 
      (particle->GetAtomicNumber()>0) ){
   if (particle->GetBaryonNumber()>0)  return true;
   else return false;
  }

   
  //  particles derived from GXIons
  if (particle->GetParticleType() == nucleus) return true;

  // proton (Hydrogen nucleus)
  if (particle->GetParticleName() == proton) return true;

  return false;
}

/////////////////
G4bool GXIonTable::IsAntiIon(const GXParticleDefinition* particle)
{
  // return true if the particle is ion

  static const G4String anti_nucleus("anti_nucleus");
  static const G4String anti_proton("anti_proton");

  // anti_neutron is not ion
  if ((particle->GetAtomicMass()>0)   && 
      (particle->GetAtomicNumber()>0) ){
   if (particle->GetBaryonNumber()<0)  return true;
   else return false;
  }

  //  particles derived from GXIons
  if (particle->GetParticleType() == anti_nucleus) return true;

  // anti_proton (Anti_Hydrogen nucleus)
  if (particle->GetParticleName() == anti_proton) return true;

  return false;
}

/////////////////
#include <algorithm>

G4bool GXIonTable::IsLightIon(const GXParticleDefinition* particle) const
{
  static const std::string names[] = { "proton", "alpha", "deuteron",
                           "triton", "He3"};

   // return true if the particle is pre-defined ion
  return std::find(names, names+5, particle->GetParticleName())!=names+5;
} 

G4bool GXIonTable::IsLightAntiIon(const GXParticleDefinition* particle) const
{
  static const std::string names[] = { "anti_proton", "anti_alpha", "anti_deuteron",
                           "anti_triton", "anti_He3"};

   // return true if the particle is pre-defined ion
  return std::find(names, names+5, particle->GetParticleName())!=names+5;
} 

/////////////////
GXParticleDefinition* GXIonTable::GetLightIon(G4int Z, G4int A) const
{
  // returns pointer to pre-defined ions
  const GXParticleDefinition* ion=0;
  if ( (Z<=2) ) {
#ifndef G4MULTITHREADED
      //In sequential use lazy-initialization
      lightions::Init();
#endif
    if ( (Z==1)&&(A==1) ) {
      ion = lightions::p_proton;
    } else if ( (Z==1)&&(A==2) ) {
        ion = lightions::p_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      ion = lightions::p_triton;
    } else if ( (Z==2)&&(A==4) ) {
      ion = lightions::p_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      ion = lightions::p_He3;
    }
  }
  return const_cast<GXParticleDefinition*>(ion);
}

/////////////////
GXParticleDefinition* GXIonTable::GetLightAntiIon(G4int Z, G4int A) const
{
  // returns pointer to pre-defined ions 
  const GXParticleDefinition* ion=0;
  if ( (Z<=2) ) {
#ifndef G4MULTITHREADED
      //In sequential use lazy-initialization
    antilightions::Init();
#endif
    if ( (Z==1)&&(A==1) ) {
      ion = antilightions::p_proton;
    } else if ( (Z==1)&&(A==2) ) {
      ion = antilightions::p_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      ion = antilightions::p_triton;
    } else if ( (Z==2)&&(A==4) ) {
      ion = antilightions::p_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      ion = antilightions::p_He3;
    }
  }
  return const_cast<GXParticleDefinition*>(ion);
}


/////////////////
// -- GetNucleusMass/GetIonMass ---
/////////////////
G4double  GXIonTable::GetNucleusMass(G4int Z, G4int A, G4int LL, G4int lvl) const
{
  if ( (A<1)  || (Z<0) || (LL<0) || (lvl<0) || (lvl>9) ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "GXIonTable::GetNucleusMass() : illegal atomic number/mass " 
             << " Z =" << Z << "  A = " << A  
	     << " L = " << LL << " lvl = " << lvl << G4endl;
    }
#endif
    G4Exception( "GXIonTable::GetNucleusMass()","PART107",
		 EventMustBeAborted, "illegal atomic number/mass");
    return -1.0;
  }
  
  G4double mass;
  if (LL == 0) {
    // calculate nucleus mass
    const GXParticleDefinition* ion=GetLightIon(Z, A);
    
    if (ion!=0) {
      mass = ion->GetPDGMass();
    } else {
      // use G4NucleiProperties::GetNuclearMass
      mass = GXNucleiProperties::GetNuclearMass(A, Z);
    }
    
    // Isomer
    if ( lvl>0 ) {
      //TK for GXBERT
      G4cout << "GXIonTable does not support GetNuclearMass(G4int Z, G4int A, G4int LL, G4int lvl) for lvl != 0 (ExcitedIsopmer)." << G4endl , abort();
      /*
      // -- loop over all particles in Ion table
      G4int encoding=GetNucleusEncoding(Z, A);
      G4IonList::iterator i = fIonList->find(encoding);
      G4bool isFound =false;
      for( ;i != fIonList->end() ; i++) {
	ion = i->second;
	if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
	// excitation level
	if ( ((const GXIons*)(ion))->GetIsomerLevel() == lvl) {
	  isFound = true;
	  break;
	}
      }
      if (isFound) {
	// return existing isomer mass
	mass = ion->GetPDGMass();
      } else {
	// Find isomer from IsotopeTable
	const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, lvl);
	if (fProperty !=0 ) mass += fProperty->GetEnergy();
      }
      */
    }

  } else {
    //TK for GXBERT
    G4cout << "GXIonTable does not support GetNuclearMass(G4int Z, G4int A, G4int LL, G4int lvl) for LL != 0 (HyperNuclei)." << G4endl , abort();
    //mass = G4HyperNucleiProperties::GetNuclearMass(A, Z, LL);
  }
  return mass;
}

/////////////////
// -- Methods for handling conatiner  ---
/////////////////

void GXIonTable::clear()
{
  if (GXParticleTable::GetParticleTable()->GetReadiness()) {
    G4Exception("GXIonTable::clear()",
		"PART116", JustWarning,
		"No effects because readyToUse is true.");    
    return;
  }

#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "GXIonTable::Clear() : number of Ion regsitered =  "; 
      G4cout << fIonList->size() <<  G4endl;
    }
#endif
  fIonList->clear();
}

void GXIonTable::Insert(const GXParticleDefinition* particle)
{
  if (!IsIon(particle)) return;
  if (Contains(particle)) return;
 
  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL); // encoding of the groud state 

  // regsiter the ion with its encoding of the groud state  
  fIonListShadow->insert( std::pair<const G4int, const GXParticleDefinition*>(encoding, particle) );

}

void GXIonTable::InsertWorker(const GXParticleDefinition* particle)
{
  if(!particle) return;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4bool found = false;
  if (encoding !=0 ) {
    G4IonList::iterator i = fIonList->find(encoding);
    for( ;i != fIonList->end() ; i++) {
      if (particle == i->second ) {
	found  = true;
	break;
      }
    }
  }
  if(found) return;
 
  // regsiter the ion with its encoding of the groud state  
  fIonList->insert( std::pair<const G4int, const GXParticleDefinition*>(encoding, particle) );

}

/////////////////
void GXIonTable::Remove(const GXParticleDefinition* particle)
{
  if(!particle) return;
  if (GXParticleTable::GetParticleTable()->GetReadiness()) {
//TK comments out for GXBERT
//    G4StateManager* pStateManager = G4StateManager::GetStateManager();
//    G4ApplicationState currentState = pStateManager->GetCurrentState();
//    if (currentState != G4State_PreInit) {
//      G4String msg = "Request of removing ";
//      msg += particle->GetParticleName();  
//      msg += " has No effects other than Pre_Init";
//      G4Exception("GXIonTable::Remove()",
//		  "PART117", JustWarning, msg);
//      return;
//    } else {
#ifdef G4VERBOSE
//      if (GetVerboseLevel()>0){
//	G4cout << particle->GetParticleName()
//	       << " will be removed from the IonTable " << G4endl;
//      }
#endif
//    }
  }

  if (IsIon(particle)) {
    G4int Z = particle->GetAtomicNumber();
    G4int A = particle->GetAtomicMass();  
    G4int LL = particle->GetQuarkContent(3);  //strangeness 
    G4int encoding=GetNucleusEncoding(Z, A, LL);
    if (encoding !=0 ) {
      G4IonList::iterator i = fIonListShadow->find(encoding);
      for( ;i != fIonListShadow->end() ; i++) {
	if (particle == i->second) {
	  fIonListShadow->erase(i);
	  break;
	}
      }
    }
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "GXIonTable::Remove :" << particle->GetParticleName() 
             << " is not ions" << G4endl; 
    }
#endif
  }
  
}



/////////////////
// -- Dump Information 
/////////////////
void GXIonTable::DumpTable(const G4String &particle_name) const
{
  const GXParticleDefinition* ion;
  G4IonList::iterator idx;
  for (idx = fIonList->begin(); idx!= fIonList->end(); ++idx) {
    ion = idx->second;
    if (( particle_name == "ALL" ) || (particle_name == "all")){
      ion->DumpTable();
    } else if ( particle_name == ion->GetParticleName() ) {
      ion->DumpTable();
    }
  }
}

/////////////////
const G4String GXIonTable::elementName[] = {
  "H",                                                                                                "He", 
  "Li", "Be",                                                             "B",  "C",  "N",  "O", "F", "Ne", 
  "Na", "Mg",                                                             "Al", "Si", "P", "S", "Cl", "Ar", 
  "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
  "Rb", "Sr", "Y", "Zr", "Nb", "Mo","Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", 
  "Cs", "Ba", 
              "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
                   "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
  "Fr", "Ra", 
              "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
              "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", 
              "Cp", "Uut", "Fl","Uup","Lv","Uus","Uuo"
};


/////////////////
G4int GXIonTable::GetVerboseLevel() const
{
  return GXParticleTable::GetParticleTable()->GetVerboseLevel();
}

/////////////////
void  GXIonTable::AddProcessManager(GXParticleDefinition* ion)
{
  // check State and do not attach process managaer in event loop
//  G4StateManager* pStateManager = G4StateManager::GetStateManager();
//  G4ApplicationState currentState = pStateManager->GetCurrentState();
//  if (currentState == G4State_EventProc) return;
//  {
//    if (n_error<10)
//    {
//      G4cout << "Defining process manager for " << ion->GetParticleName() << G4endl;
//      G4Exception("GXIonTable::AddProcessManager()", "PART130", JustWarning,
//	"Defining process manager during an event loop is thread unsafe and will be dropped from the next release.");
//      n_error +=1;
//    }
//    return;
//  }

  if(ion->IsGeneralIon()) {

    // check whether GenericIon has processes
    GXParticleDefinition* genericIon = 
      GXParticleTable::GetParticleTable()->GetGenericIon();

    G4ProcessManager* pman=0;
    if (genericIon!=0) pman = genericIon->GetProcessManager();
    if ((genericIon ==0) || (genericIon->GetParticleDefinitionID() < 0) || (pman==0)){
      G4cout << "GXIonTable::AddProcessManager() : can not create ion of  " 
             << ion->GetParticleName()
             << "  because GenericIon is not available!!" <<   G4endl;
      G4Exception( "GXIonTable::AddProcessManager()","PART105", FatalException, 
                   "Can not create ions because GenericIon is not available");
      return;
    }
  
    ////////  ion->SetProcessManager(pman);
    //TK for GXBERT
    //ion->SetParticleDefinitionID(genericIon->GetParticleDefinitionID());
  }
  else {

    // is this a MuonicAtom ?

G4cout << "is this a MuonicAtom ?" << G4endl, abort();
/*
    G4MuonicAtom* muatom = dynamic_cast<G4MuonicAtom*> (ion);

    if ( muatom != nullptr ) {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout << "GXIonTable::AddProcessManager() : MuonicAtom dynamic_cast succeeded for " 
               << ion->GetParticleName() << G4endl;
      }
#endif
      // check whether GenericMuonicAtom has processes
      GXParticleDefinition* genericMA = 
        GXParticleTable::GetParticleTable()->GetGenericMuonicAtom();

      G4ProcessManager* pman = nullptr;
      if (genericMA != nullptr) pman = genericMA->GetProcessManager();
      if ((genericMA == nullptr) || (genericMA->GetParticleDefinitionID() < 0) || (pman==nullptr)){
        G4cout << "GXIonTable::AddProcessManager() : can not create MuonicAtom  " 
               << ion->GetParticleName()
               << "  because GenericMuonicAtom is not available!!" <<   G4endl;
        G4Exception( "GXIonTable::AddProcessManager()","PART106", FatalException, 
                     "Can not create MuonicAtoms because GenericMuonicAtom is not available");
        return;
      }
  
      ////////  ion->SetProcessManager(pman);
      ion->SetParticleDefinitionID(genericMA->GetParticleDefinitionID());
      
    }
    else {
      G4cout << "GXIonTable::AddProcessManager() : can not create  " 
             << ion->GetParticleName()
             << "  because of unsupported particle type !!" <<   G4endl;
      G4Exception( "GXIonTable::AddProcessManager()","PART107", FatalException, 
                   "Can not create particle");
      return;
    }
*/

  }
  return;
}

#include <vector>     

////////////////////
GXParticleDefinition* GXIonTable::GetParticle(G4int index) const
{
  if ( (index >=0) && (index < Entries()) ) {
    G4IonList::iterator idx = fIonList->begin();
    G4int counter = 0;
    while( idx != fIonList->end() ){// Loop checking, 09.08.2015, K.Kurashige
      if ( counter == index ) {
	return const_cast<GXParticleDefinition*>(idx->second);
      }
      counter++;
      idx++;
    }
  } 
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1){
    G4cout << " GXIonTable::GetParticle"
           << " invalid index (=" << index << ")" 
	   << " entries = " << Entries() << G4endl;
  }
#endif
  return 0;
}

////////////////////
G4bool  GXIonTable::Contains(const GXParticleDefinition* particle) const
{
  if (!IsIon(particle)) return false;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4bool found = false;
  if (encoding !=0 ) {
    G4IonList::iterator i = fIonListShadow->find(encoding);
    for( ;i != fIonListShadow->end() ; i++) {
      if (particle == i->second ) {
	found  = true;
	break;
      }
    }
  }
  return found;
}

////////////////////
G4int GXIonTable::Entries() const
{
  return fIonList->size();
}

////////////////////
G4int GXIonTable::size() const
{
  return fIonList->size();
}


////////////////////
GXParticleDefinition* GXIonTable::FindIonInMaster(G4int Z, G4int A, G4double E, 
                            GXIons::GXFloatLevelBase flb, G4int /*J*/)
{
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const GXParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    // excitation level
    G4double anExcitaionEnergy = ((const GXIons*)(ion))->GetExcitationEnergy();
    //TK for GXBERT
    //if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
    if (std::fabs(E - anExcitaionEnergy) < 1.0*CLHEP::eV ) {
      if(((const GXIons*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound ){ 
    return const_cast<GXParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
GXParticleDefinition* GXIonTable::FindIonInMaster(G4int Z, G4int A, G4int LL, G4double E,
                         GXIons::GXFloatLevelBase flb, G4int J)
{
  if (LL==0) return FindIon(Z,A,E,flb,J);
  
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const GXParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL, 0.0, 0);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if(  ion->GetQuarkContent(3) != LL) break;
    // excitation level
    G4double anExcitaionEnergy = ((const GXIons*)(ion))->GetExcitationEnergy();
    //TK for GXBERT
    //if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
    if (std::fabs(E - anExcitaionEnergy) < 1.0*CLHEP::eV ) {
      if(((const GXIons*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound ){ 
    return const_cast<GXParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
GXParticleDefinition* GXIonTable::FindIonInMaster(G4int Z, G4int A, G4int lvl)
{
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const GXParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    // excitation level
    if ( ((const GXIons*)(ion))->GetIsomerLevel() == lvl) {
      isFound = true;
      break;
    } 
  }

  if ( isFound ){ 
    return const_cast<GXParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
GXParticleDefinition* GXIonTable::FindIonInMaster(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return FindIon(Z,A,lvl);
  
  // Search ions with A, Z ,E, lvl
  const GXParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if ( ion->GetQuarkContent(3) != LL) break;
    // excitation level
    if ( ((const GXIons*)(ion))->GetIsomerLevel() == lvl) {
      isFound = true;
      break;
    }
  }

  if ( isFound ){ 
    return const_cast<GXParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4double GXIonTable::GetLifeTime(const GXParticleDefinition* particle) const
{
  //TK for GXBERT 
  if ( !(particle->IsGeneralIon()) ) {
     return particle->GetPDGLifeTime();
  } else {
     G4cout << "GXIonTable::GetLifeTime(const GXParticleDefinition* particle) does not supported for GeneralIon" << G4endl; abort();
  }
/*
  //if(!(particle->IsGeneralIon())) return particle->GetPDGLifeTime();

  //const GXIons* ion = static_cast<const GXIons*>(particle);
  //G4int Z = ion->GetAtomicNumber();
  //G4int A = ion->GetAtomicMass();
  //G4double E = ion->GetExcitationEnergy();

  if((particle->IsGeneralIon()) && !pNuclideTable)
  {
   G4Exception("GXIonTable::GetLifeTime()","ParticleIon1001",FatalException,
               "Method is invoked before GXIonTable is initialized.");
   //return 0.;
  } //else {
   //G4IsotopeProperty* isoP = pNuclideTable->GetIsotope(Z,A,E);
    //if(!isoP) return -1001.0;
    //return isoP->GetLifeTime();
  //}
  return particle->GetPDGLifeTime();
*/
}

////////////////////
G4double GXIonTable::GetLifeTime(G4int Z, G4int A, G4double E, char flbChar) const
{ return GetLifeTime(Z,A,E,GXIons::FloatLevelBase(flbChar)); }

////////////////////
G4double GXIonTable::GetLifeTime(G4int Z, G4int A, G4double E,
             GXIons::GXFloatLevelBase flb) const
{
  //TK for GXBERT 
  G4cout << "GXIonTable::GetLifeTime(G4int Z, G4int A, G4double E, GXIons::GXFloatLevelBase flb) is not supported" << G4endl; abort();
  return -1001.0;
/*
  G4double life = -1001.0;
  const G4IsotopeProperty* fProperty = FindIsotope(Z, A, E, flb);
  if( fProperty !=0 ) life = fProperty->GetLifeTime();
  return life;
*/
}
