
#include "VecCore/VecCore"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzConvertor.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"
#include "GXLorentzConvertor.hh"
#include "GXProton.hh"
#include "GXNeutron.hh"
#include "LorentzVector.hh"

using namespace gxbert;

int main(int argc, char* argv[]) {

  int debugLevel = 5;

  //.. default run with a fixed target and energy
  double energy =  1.3E+3; // in MeV
  int nEvents   = 1; //1024*64;
  if(argc >= 2) energy  = atof(argv[1]);
  if(argc >= 3) nEvents = atoi(argv[2]);

  //.. timer
  // Timer<nanoseconds> timer;
  // double timeElapsed = 0.;

  //.. prepare bullets
  GXTrackHandler *bulletHandler = new GXTrackHandler(nEvents);
  bulletHandler->GenerateRandomTracks(nEvents, energy, energy);
  GXTrack* aosBullets = bulletHandler->GetAoSTracks();
  //printf("Number of tracks = %d\n", soaBullets.size);

  //.. prepare targets
  GXTrackHandler *targetHandler = new GXTrackHandler(nEvents);
  targetHandler->GenerateRandomTracks(nEvents, energy, energy);
  GXTrack* aosTargets = targetHandler->GetAoSTracks();

  //==== Reference: G4LorentzConvertor functionalities
  double g4sumEscm = 0.;
  double g4sumEkin = 0.;
  double g4sumP2 = 0.;
  G4LorentzConvertor g4conv;
  g4conv.setVerbose(debugLevel);
  for(size_t i = 0; i < nEvents; ++i) {
    const GXTrack &trk1 = aosBullets[i];
    G4LorentzVector lorvec1(trk1.px, trk1.py, trk1.pz, trk1.E);
    G4InuclElementaryParticle* particle1 = new G4InuclElementaryParticle(lorvec1, GXProton::Definition());

    const GXTrack &trk2 = aosTargets[i];
    G4LorentzVector lorvec2(trk2.px, trk2.py, trk2.pz, trk2.E);
    G4InuclElementaryParticle* particle2 = new G4InuclElementaryParticle(lorvec2, GXNeutron::Definition());
    g4conv.setBullet(particle1);
    g4conv.setTarget(particle2);

    g4conv.toTheCenterOfMass();
    g4sumEscm += g4conv.getTotalSCMEnergy();
    g4sumEkin += g4conv.getKinEnergyInTheTRS();

    G4LorentzVector plab = g4conv.backToTheLab(particle1->getMomentum());
    g4sumP2 += plab.mag2();
  }

  std::cout<<"GXBert results: "
	   <<"  sumEscm = "<< g4sumEscm
	   <<"\t  sumEkin = "<< g4sumEkin
	   <<"\t  sumP2 = "<< g4sumP2
	   <<"\n";


  //=== Templated converter in scalar mode
  double gsSumEscm = 0.;
  double gsSumEkin = 0.;
  double gsSumP2 = 0.;
  GXLorentzConvertor<double> gsconv;
  gsconv.setVerbose(debugLevel);

  GXTrack_v &soaBullets = bulletHandler->GetSoATracks();
  GXTrack_v &soaTargets = targetHandler->GetSoATracks();
  LorentzVector<double> fourvec;
  for(size_t i = 0; i < nEvents; i++) {
    soaBullets.getFourMomentum(i, fourvec);
    gsconv.setBullet(fourvec);

    soaTargets.getFourMomentum(i, fourvec);
    gsconv.setTarget(fourvec);

    gsconv.toTheCenterOfMass();
    gsSumEscm += g4conv.getTotalSCMEnergy();
    gsSumEkin += g4conv.getKinEnergyInTheTRS();

    // LorentzVector<double> plab = g4conv.backToTheLab(particle1->getMomentum());
    // gsSumP2 += plab.mag2();    
  }

  std::cout<<"Scalar results: "
	   <<"  sumEscm = "<< gsSumEscm
	   <<"\t  sumEkin = "<< gsSumEkin
	   <<"\t  sumP2 = "<< gsSumP2
	   <<"\n";


  //=== Vectorized converter
  //.. vector size 
  // using Double_v = typename vecCore::backend::VcVector::Double_v;
  // int vsize = vecCore::VectorSize<Double_v>();

  // double gxSumEscm = 0.;
  // double gxSumEkin = 0.;
  // double gxSumP2 = 0.;
  // GXLorentzConvertor<Double_v> gxconv;
  // GXTrack_v const& soaBullets = bulletHandler->GetSoATracks();
  // GXTrack_v const& soaTargets = targetHandler->GetSoATracks();
  // for(size_t i = 0; i < nEvents; i+=vsize) {
  //   //
  // }

  //.. cleanup
  delete bulletHandler;
  delete targetHandler;
}
