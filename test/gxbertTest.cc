#include "GXCascadeInterface.hh"
#include "GXDynamicParticle.hh"

#include "GXParticleTable.hh"
#include "GXProton.hh"
#include "GXNeutron.hh"
#include "GXPionPlus.hh"
#include "GXPionMinus.hh"
#include "GXPionZero.hh"
#include "GXGamma.hh"
#include "GXKaonPlus.hh"
#include "GXKaonMinus.hh"
#include "GXKaonZero.hh"
#include "GXKaonZeroLong.hh"
#include "GXKaonZeroShort.hh"
#include "GXAntiKaonZero.hh"
#include "GXLambda.hh"
#include "GXSigmaPlus.hh"
#include "GXSigmaZero.hh"
#include "GXSigmaMinus.hh"
#include "GXXiZero.hh"
#include "GXXiMinus.hh"
#include "GXOmegaMinus.hh"
#include "GXDeuteron.hh"
#include "GXTriton.hh"
#include "GXHe3.hh"
#include "GXAlpha.hh"

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Random/Randomize.h>

#include<map>
#include <iostream>
#include <numeric>

void output( G4int , std::vector< std::pair<G4double,G4double> >& );

int main(int argc , char* const argv[] ) {

   std::string s_proj = argv[1];
   std::string s_targ = argv[2];
   G4double energy = atof(argv[3]);
   int event_per_job(atoi(argv[4]));

   std::map< std::string , GXNucleus* > mNucleus;
   mNucleus.insert( std::pair< std::string , GXNucleus* > ("C", new GXNucleus(   12 ,  6) ) );
   mNucleus.insert( std::pair< std::string , GXNucleus* > ("Al", new GXNucleus(  27 , 13) ) );
   mNucleus.insert( std::pair< std::string , GXNucleus* > ("Fe", new GXNucleus(  56 , 26) ) );
   mNucleus.insert( std::pair< std::string , GXNucleus* > ("In", new GXNucleus( 115 , 49) ) );
   mNucleus.insert( std::pair< std::string , GXNucleus* > ("Pb", new GXNucleus( 208 , 82) ) );

   std::vector<GXParticleDefinition*> v_PD;
   v_PD.push_back( GXProton::Proton() );
   v_PD.push_back( GXNeutron::Neutron() );
   v_PD.push_back( GXGamma::Gamma() );
   v_PD.push_back( GXKaonMinus::KaonMinus() );
   v_PD.push_back( GXKaonPlus::KaonPlus() );
   v_PD.push_back( GXLambda::Lambda() );
   v_PD.push_back( GXPionMinus::PionMinus() );
   v_PD.push_back( GXPionPlus::PionPlus() );
   v_PD.push_back( GXSigmaMinus::SigmaMinus() );
   v_PD.push_back( GXSigmaPlus::SigmaPlus() );
   v_PD.push_back( GXOmegaMinus::OmegaMinus() );
   v_PD.push_back( GXXiMinus::XiMinus() );


   GXCascadeInterface* cascadeInterface = new GXCascadeInterface;
   cascadeInterface->SetVerboseLevel(0);
   GXHadFinalState* result = new GXHadFinalState;

   GXParticleDefinition* pd = GXParticleTable::GetParticleTable()->FindParticle( s_proj );

   GXNucleus nucleus = *(mNucleus.find( s_targ )->second);

   //prepair projectile
   //GXDynamicParticle* dp = new GXDynamicParticle( pd , G4ThreeVector(0,0,1) , energy[ie]*CLHEP::GeV );
   //GXHadProjectile projectile(*dp);
   GXHadProjectile projectile( pd , energy*CLHEP::GeV );
   std::cerr<<" gxbert bullet: 4mom=("<< projectile.fourMomentum.x() <<"; "<< projectile.fourMomentum.y() <<"; "<< projectile.fourMomentum.z() <<"; "<< projectile.fourMomentum.t() <<") and mass="<< projectile.fourMomentum.mag() <<"\n";

   //prepair container for secondary neutrons
   //                           ke      z:mu
   std::vector< std::pair<G4double,G4double> > neutrons;

   for ( G4int ij = 0 ; ij != event_per_job; ij++ ) { //Loop events(reactions)
   //make a reaction 
      result = cascadeInterface->ApplyYourself( projectile , nucleus );

      //for ( G4int i = 0 ; i != result->GetNumberOfSecondaries() ; i++ ) { //Loop secondary particles 
      for ( G4int i = 0 ; i != (G4int)result->GetSecondaries().size() ; i++ ) { //Loop secondary particles 
         //GXDynamicParticle* dp = result->GetSecondary(i)->GetParticle();
         GXDynamicParticle* dp = result->GetSecondaries().at(i);
         //register secondary neutrons 
         if (  dp->GetDefinition() == GXNeutron::Neutron() ) neutrons.push_back( std::pair<G4double,G4double>( dp->GetKineticEnergy() , dp->GetMomentumDirection().z() ) );

         //if ( cascadeInterface->GetVerboseLevel() > 0 ) {
/*
                     G4cout << "RESULT"
                     << " " << ip << " " << it << " " << ie << " " << ij << " " << i
                     << " " << dp->GetDefinition()->GetParticleName()
                     << " " << dp->GetKineticEnergy()
                     << " " << dp->Get4Momentum()
                     //<< " " << dp->GetMomentumDirection()
                            << G4endl;
*/
         //}
      } //Loop seconary particles 
  }//Loop events(reactions)

  output( event_per_job , neutrons );

}

void output( G4int NbOfEvents , std::vector< std::pair<G4double,G4double> >& neutrons ) {

   //Re use codes from SATIF_CYL benchmark
   G4double pi=CLHEP::pi;

   G4int iemin=-110; //1E-5eV
   G4int iemax=100;  //10,000TeV 
   std::map<G4int,G4int> mAngNumber;
   for ( G4int i = 0 ; i <= 180 ; i++ ) {
      mAngNumber.insert( std::pair<G4int,G4int>(i,0) );
   }

   std::map< G4int,std::map<G4int,G4int> > mAngEneNumber;
   for ( G4int iang = 0 ; iang <= 180 ; iang++ ) {
      std::map<G4int,G4int> mEneNumber;
      for ( G4int ie = iemin ; ie <= iemax ; ie++ ) {
         mEneNumber.insert( std::pair<G4int,G4int>(ie,0) );
      }
      mAngEneNumber.insert( std::pair<G4int,std::map<G4int,G4int> >(iang,mEneNumber) );
   }

   for ( std::vector< std::pair<G4double,G4double> >::iterator
         it = neutrons.begin() ; it != neutrons.end() ; it++ ) {
      G4int i = floor ( acos(it->second)*180/pi + 0.5 );
      if ( it->first > 20*CLHEP::MeV ) mAngNumber.find( i )->second++;

      G4double log10_e = log10( it->first/CLHEP::MeV );
      G4int ie = floor( log10_e*10 );
      mAngEneNumber.find( i )->second.find( ie )->second++;
   }

   G4double dsr0=2*pi*(cos(0/180*pi)-cos(0.5/180*pi));
   G4double dsr15=2*pi*(cos(14.5/180*pi)-cos(15.5/180*pi));
   G4double dsr30=2*pi*(cos(29.5/180*pi)-cos(30.5/180*pi));
   G4double dsr45=2*pi*(cos(44.5/180*pi)-cos(45.5/180*pi));
   G4double dsr60=2*pi*(cos(59.5/180*pi)-cos(60.5/180*pi));
   G4double dsr90=2*pi*(cos(89.5/180*pi)-cos(90.5/180*pi));
   G4double dsr120=2*pi*(cos(119.5/180*pi)-cos(120.5/180*pi));
   G4double dsr150=2*pi*(cos(149.5/180*pi)-cos(150.5/180*pi));

   G4int ie_begin = 0;
   G4int ie_end = 60;

   //double differentional
   G4cout << "                lower_side_energy 0deg 15deg 30deg 45deg 60deg 90deg 120deg 150deg AngInt"
          << G4endl;
   for ( G4int ie = ie_begin ; ie <= ie_end ; ie++ ) {

      G4double e = pow( 10.0 , (1.0*ie) / 10 );
      G4double de = pow( 10.0 , (1.0*ie+1) / 10 ) - e;
      G4double factor = 1.0/de/NbOfEvents;

      G4int isum = 0;
      for ( G4int iang = 0 ; iang <= 180 ; iang++ ) isum += mAngEneNumber.find(iang)->second.find(ie)->second;

      G4cout << "RESULT:SATIF_CYL.rev1"
      << " " << e
      << " " << mAngEneNumber.find(0)->second.find(ie)->second*factor/dsr0
      << " " << mAngEneNumber.find(15)->second.find(ie)->second*factor/dsr15
      << " " << mAngEneNumber.find(30)->second.find(ie)->second*factor/dsr30
      << " " << mAngEneNumber.find(45)->second.find(ie)->second*factor/dsr45
      << " " << mAngEneNumber.find(60)->second.find(ie)->second*factor/dsr60
      << " " << mAngEneNumber.find(90)->second.find(ie)->second*factor/dsr90
      << " " << mAngEneNumber.find(120)->second.find(ie)->second*factor/dsr120
      << " " << mAngEneNumber.find(150)->second.find(ie)->second*factor/dsr150
      << " " << isum*factor
             << G4endl;

   }

   //direction differentional [/dOmega
   G4cout << "                EergyIntFluxN>20MeV at 0deg 15deg 30deg 45deg 60deg 90deg 120deg 150deg TOT " << G4endl;
   G4cout << "RESULT:SATIF_CYL.rev1"
   << " IntE>20MeV"
   << " " << mAngNumber.find(   0 )->second/dsr0/NbOfEvents
   << " " << mAngNumber.find(  15 )->second/dsr15/NbOfEvents
   << " " << mAngNumber.find(  30 )->second/dsr30/NbOfEvents
   << " " << mAngNumber.find(  45 )->second/dsr45/NbOfEvents
   << " " << mAngNumber.find(  60 )->second/dsr60/NbOfEvents
   << " " << mAngNumber.find(  90 )->second/dsr90/NbOfEvents
   << " " << mAngNumber.find( 120 )->second/dsr120/NbOfEvents
   << " " << mAngNumber.find( 150 )->second/dsr150/NbOfEvents
   << " " << std::accumulate( std::begin( mAngNumber )
                            , std::end( mAngNumber )
                            , 0
                            , [] ( G4int value , const std::map<G4int,G4int>::value_type& p){ return value +p.second;} ) /1.0/NbOfEvents
          << G4endl;
}

