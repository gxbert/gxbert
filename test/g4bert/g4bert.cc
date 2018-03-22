#include "G4CascadeInterface.hh"
#include "G4DynamicParticle.hh"
#include "G4ProcessManager.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ParticleTable.hh"

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Random/Randomize.h>

#include <numeric>
#include <iostream>

void output( G4int , std::vector< std::pair<G4double,G4double> >& );

int main( int argc , char* const argv[] ) {

   std::string s_proj = argv[1];
   std::string s_targ = argv[2];
   std::string s_energy = argv[3];

   int event_per_job(1024);
   //int event_per_job(1024*64);

   //prepair energy 
   std::map< std::string , double > mEnergy;
   mEnergy.insert( std::pair< std::string , double > ( "800MeV" , 0.8 ) );
   mEnergy.insert( std::pair< std::string , double > ( "1.5GeV" , 1.5 ) );
   mEnergy.insert( std::pair< std::string , double > ( "3GeV"   , 3.0 ) );

   //prepair target 
   std::map< std::string , G4Nucleus* > mNucleus;
   mNucleus.insert( std::pair< std::string , G4Nucleus* > ("C", new G4Nucleus(  12 ,   6 ) ) );
   mNucleus.insert( std::pair< std::string , G4Nucleus* > ("Al", new G4Nucleus( 27 ,  13 ) ) );
   mNucleus.insert( std::pair< std::string , G4Nucleus* > ("Fe", new G4Nucleus( 56 ,  26 ) ) );
   mNucleus.insert( std::pair< std::string , G4Nucleus* > ("In", new G4Nucleus( 115 , 49 ) ) );
   mNucleus.insert( std::pair< std::string , G4Nucleus* > ("Pb", new G4Nucleus( 209 , 82 ) ) );

   //prepair particle definition of projectile 
   std::vector<G4ParticleDefinition*> v_PD;
   v_PD.push_back( G4Gamma::Gamma() );
   v_PD.push_back( G4Proton::Proton() );
   v_PD.push_back( G4Neutron::Neutron() );
   v_PD.push_back( G4KaonMinus::KaonMinus() );
   v_PD.push_back( G4KaonPlus::KaonPlus() );
   v_PD.push_back( G4Lambda::Lambda() );
   v_PD.push_back( G4PionMinus::PionMinus() );
   v_PD.push_back( G4PionPlus::PionPlus() );
   v_PD.push_back( G4SigmaMinus::SigmaMinus() );
   v_PD.push_back( G4SigmaPlus::SigmaPlus() );
   v_PD.push_back( G4OmegaMinus::OmegaMinus() );
   v_PD.push_back( G4XiMinus::XiMinus() );
   
   //prepair particles in case they are needed in reaction
   G4BosonConstructor  pBosonConstructor;
   pBosonConstructor.ConstructParticle();
   G4LeptonConstructor pLeptonConstructor;
   pLeptonConstructor.ConstructParticle();
   G4MesonConstructor pMesonConstructor;
   pMesonConstructor.ConstructParticle();
   G4BaryonConstructor pBaryonConstructor;
   pBaryonConstructor.ConstructParticle();
   G4IonConstructor pIonConstructor;
   pIonConstructor.ConstructParticle();

   //register dummy process manager to avoid error
   G4GenericIon::GenericIon()->SetProcessManager(new G4ProcessManager( G4GenericIon::GenericIon() ));

   //prepair BERT cascade 
   G4CascadeInterface* cascadeInterface = new G4CascadeInterface;
   cascadeInterface->SetVerboseLevel(0);

   G4HadFinalState* result = new G4HadFinalState;
   //G4DynamicParticle* dp = new G4DynamicParticle( G4Proton::Proton() , G4ThreeVector(0,0,1) , 1*CLHEP::GeV );
   //G4HadProjectile projectile(*dp);
   //G4Nucleus nucleus(12,6);

   G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->FindParticle( s_proj );

   G4Nucleus nucleus = *(mNucleus.find( s_targ )->second);

   G4double energy = mEnergy.find( s_energy )->second;

   //prepair projectile
   G4DynamicParticle* dp = new G4DynamicParticle( pd , G4ThreeVector(0,0,1) , energy*CLHEP::GeV );
   G4HadProjectile projectile(*dp);

   //prepair container for secondary neutrons
   //                           ke      z:mu
   std::vector< std::pair<G4double,G4double> > neutrons;
           
   for ( G4int ij = 0 ; ij != event_per_job; ij++ ) { //Loop events(reactions)
      //make a reaction 
      result = cascadeInterface->ApplyYourself( projectile , nucleus );

      for ( G4int i = 0 ; i != result->GetNumberOfSecondaries() ; i++ ) { //Loop secondary particles 
         G4DynamicParticle* dp = result->GetSecondary(i)->GetParticle();
         //register secondary neutrons 
         if (  dp->GetDefinition() == G4Neutron::Neutron() ) neutrons.push_back( std::pair<G4double,G4double>( dp->GetKineticEnergy() , dp->GetMomentumDirection().z() ) );

          if ( cascadeInterface->GetVerboseLevel() > 0 ) {
             G4cout << "RESULT" 
             << " " << dp->GetDefinition()->GetParticleName()
             << " " << dp->GetKineticEnergy()
             << " " << dp->Get4Momentum()
             << " " << dp->GetMomentumDirection()
                    << G4endl;
           }
       } //Loop seconary particles 
   }//Loop events(reactions)

   output( event_per_job , neutrons );

   delete dp;
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
