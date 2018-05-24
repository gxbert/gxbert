
#include "VecCore/VecCore"
#include "VecCore/Timer.h"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclCollider.hh"
#include "G4CollisionOutput.hh"
#include "LorentzVector.hh"
#include "GXParticleDefinition.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"
//#include "GXInuclCollider.hh"

#include "GXProton.hh"
#include "GXNeutron.hh"
#include "G4InuclParticleNames.hh"

#include <iostream>

using namespace gxbert;
using std::cerr;

//.. default globals
constexpr double pmass = 0.938272013;  // proton mass in GeV
int nReps        = 1;
int nEvents      = 1; // 1024*64;
double kinEnergy = 1.500; // in GeV

int debugLevel   = 1;

bool isOutputInvalid(G4InuclParticle const &bullet, G4CollisionOutput const &output)
{
  // Quantities necessary for conditional block below
  G4int npart = output.numberOfOutgoingParticles();
  G4int nfrag = output.numberOfOutgoingNuclei();

  const GXParticleDefinition* firstOut = (npart == 0) ? 0 :
    output.getOutgoingParticles().begin()->getDefinition();

  /*
#ifdef G4CASCADE_DEBUG_INTERFACE
  // Report on all retry conditions, in order of return logic
  G4cout << " retryInelasticNucleus: numberOfTries "
	 << ((numberOfTries < maximumTries) ? "RETRY (t)" : "EXIT (f)")
	 << "\n retryInelasticNucleus: AND outputParticles "
	 << ((npart != 0) ? "NON-ZERO (t)" : "EMPTY (f)")
#ifdef G4CASCADE_COULOMB_DEV
	 << "\n retryInelasticNucleus: AND coulombBarrier (COULOMB_DEV) "
	 << (coulombBarrierViolation() ? "VIOLATED (t)" : "PASSED (f)")
	 << "\n retryInelasticNucleus: AND collision type (COULOMB_DEV) "
	 << ((npart+nfrag > 2) ? "INELASTIC (t)" : "ELASTIC (f)")
#else
	 << "\n retryInelasticNucleus: AND collsion type "
	 << ((npart+nfrag < 3) ? "ELASTIC (t)" : "INELASTIC (f)")
	 << "\n retryInelasticNucleus: AND Leading particle bullet "
	 << ((firstOut == bullet->getDefinition()) ? "YES (t)" : "NO (f)")
#endif
	 << "\n retryInelasticNucleus: OR conservation "
	 << (!balance->okay() ? "FAILED (t)" : "PASSED (f)")
	 << G4endl;
#endif
  */

  return ( ((npart != 0) &&
	    //#ifdef G4CASCADE_COULOMB_DEV
	    //(coulombBarrierViolation() && npart+nfrag > 2)
	    //#else
	      (npart+nfrag < 3 && firstOut == bullet.getDefinition())
	    //#endif
             )
	     //#ifndef G4CASCADE_SKIP_ECONS
	     //|| (!balance->okay())
	     //#endif
	 );
}

//==== Reference: G4InuclCollider functionalities
void RunG4InuclColliderTimer(GXTrack_v const& soaBullets, GXTrack_v const& soaTargets, G4CollisionOutput &output)
{
  double g4sumEscm = 0.;
  double g4sumEkin = 0.;
  double g4sumP2 = 0.;

  G4InuclCollider collider;
  collider.setVerboseLevel(debugLevel);

  G4InuclElementaryParticle bullet;
  G4InuclNuclei target;
  G4LorentzVector lorvec;
  Timer<nanoseconds> timer;
  timer.Start();
  int numberOfTries = 0;
  for(size_t irep = 0; irep < nReps; ++irep) {
    g4sumEscm = g4sumEkin = g4sumP2 = 0.;
    for(size_t i = 0; i < nEvents; ++i) {
      soaBullets.getFourMomentum(i, lorvec);
      bullet.fill(lorvec, G4InuclParticleNames::proton);

      soaTargets.getFourMomentum(i, lorvec);
      target.fill(lorvec, 208, 82); // lead: A=209, Z=82

      numberOfTries = 0;
      do {   			// we try to create inelastic interaction
	if (debugLevel > 1) cerr<<"=== Generating cascade attempt " << numberOfTries <<"\n";

	output.reset();
	collider.collide(&bullet, &target, output);

	numberOfTries++;
      } while ( numberOfTries <= 20 && isOutputInvalid(bullet, output) );
      cerr<<"***** Event "<< i <<" Ntries="<< numberOfTries <<" Final state: "<< output.numberOfOutgoingParticles() <<" hadrons + "<< output.numberOfOutgoingNuclei() <<" nuclei\n";
    }
  }
  double g4Elapsed = timer.Elapsed();

  std::cerr<<"GXBert results: "
	   <<"  sumEscm = "<< g4sumEscm
	   <<"\t  sumEkin = "<< g4sumEkin
	   <<"\t  sumP2 = "<< g4sumP2
	   <<"\t CPUtime = "<< g4Elapsed / nEvents / nReps
	   <<"\n";
}

/*
template <typename Real_v>
void RunGXInuclColliderTimer(const char* testname, GXTrack_v const& soaBullets, GXTrack_v const& soaTargets)
{
  int vsize = vecCore::VectorSize<Real_v>();
  Real_v sumEscm = 0.;
  Real_v sumEkin = 0.;
  Real_v sumP2 = 0.;
  GXInuclCollider<Real_v> conv;
  conv.setVerbose(debugLevel);

  LorentzVector<Real_v> fourvec;
  Timer<nanoseconds> timer;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sumEscm = sumEkin = sumP2 = 0.;
    for(size_t i = 0; i < nEvents; i += vsize) {
      soaTargets.getFourMomentum(i, fourvec);
      conv.setTarget(fourvec);

      soaBullets.getFourMomentum(i, fourvec);
      conv.setBullet(fourvec);

      conv.toTheCenterOfMass();
      sumEscm += conv.getTotalSCMEnergy();
      sumEkin += conv.getKinEnergyInTheTRS();

      LorentzVector<Real_v> plab = conv.backToTheLab(fourvec);
      sumP2 += plab.Mag2();
    }
  }
  double elapsed = timer.Elapsed();

  std::cerr<< testname <<" results: "
	   <<   "sumEscm = "<< vecCore::ReduceAdd(sumEscm)
	   <<"\t sumEkin = "<< vecCore::ReduceAdd(sumEkin)
	   <<"\t sumP2 = "  << vecCore::ReduceAdd(sumP2)
	   <<"\t CPUtime = "<< elapsed / nEvents / nReps
	   <<"\n";
}
*/

int main(int argc, char* argv[]) {

  if (argc >= 2) debugLevel = atoi(argv[1]);
  if (argc >= 3) nEvents = atoi(argv[2]);
  if (argc >= 4) nReps = atoi(argv[3]);
  if (argc >= 5) kinEnergy  = atof(argv[4]); // energy in GeV

  const double etot = kinEnergy + pmass;
  double pmom = sqrt(etot * etot - pmass * pmass);
  cerr<<"#evts="<< nEvents <<" x #reps="<< nReps <<" at Ek = "<< kinEnergy <<"GeV"<<" - Etot="<< etot <<"GeV and pmom="<< pmom <<"GeV\n";

  //.. timer
  Timer<nanoseconds> timer;

  //.. prepare bullets - note that the handler uses units in MeV! (since the proton mass is hardcoded)
  GXTrackHandler *bulletHandler = new GXTrackHandler(nEvents); 
  bulletHandler->SetMass( 0.938272013 );  //proton mass in GeV
  double posdir[6] = {0., 0., 0., 0., 0., 1.};
  bulletHandler->GenerateTracksAlongSameDirection(nEvents, posdir, pmom, pmom);
  GXTrack_v &soaBullets = bulletHandler->GetSoATracks();
  printf("Number of tracks = %d at Ekin=%f GeV and pMom=%f GeV\n", soaBullets.size, kinEnergy, pmom);

  //.. prepare targets
  GXTrackHandler *targetHandler = new GXTrackHandler(nEvents);
  targetHandler->SetMass( 207.976636 );  // proton mass in GeV
  targetHandler->GenerateTracksAlongSameDirection(nEvents, posdir, 0., 0.);
  GXTrack_v &soaTargets = targetHandler->GetSoATracks();
  double nucleusMass = G4InuclNuclei(0.,208,82,0.).getNucleiMass();
  for(int i=0; i<nEvents; ++i) {
    soaTargets.E[i] = nucleusMass;
  }
  
  // Using function calls for benchmarks
  G4CollisionOutput output;
  RunG4InuclColliderTimer(soaBullets, soaTargets, output);

  // benchmarks templated on types
  //RunGXInuclColliderTimer<double>("double", soaBullets, soaTargets);
  //RunGXInuclColliderTimer<vecCore::backend::VcVector::Double_v>("Double_v", soaBullets, soaTargets);

  /*
  //=== Alternate way to load arrays for vectorized converter
  Double_v px, py, pz, m, kinE;
  Real_v glSumEscm = 0.;
  Real_v glSumEkin = 0.;
  Real_v glSumP2 = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    glSumEscm = glSumEkin = glSumP2 = 0.;
    for (int i = 0; i < nEvents; i += vsize) {
      Load(px, &soaTargets.px[i]);
      Load(py, &soaTargets.py[i]);
      Load(pz, &soaTargets.pz[i]);
      Load( m, &soaTargets.m[i]);
      Load(kinE, &soaTargets.E[i]);
      simd4vec.Set(px, py, pz, kinE+m);
      gxconv.setTarget(simd4vec);//, simd4vec.Mag());

      Load(px, &soaBullets.px[i]);
      Load(py, &soaBullets.py[i]);
      Load(pz, &soaBullets.pz[i]);
      Load( m, &soaBullets.m[i]);
      Load(kinE, &soaBullets.E[i]);
      simd4vec.Set(px, py, pz, kinE+m);
      gxconv.setBullet(simd4vec);//, simd4vec.Mag());

      gxconv.toTheCenterOfMass();
      glSumEscm += gxconv.getTotalSCMEnergy();
      glSumEkin += gxconv.getKinEnergyInTheTRS();

      LorentzVector<Real_v> plab = gxconv.backToTheLab(simd4vec);
      glSumP2 += plab.Mag2();
    }
  }
  double glElapsed = timer.Elapsed();

  std::cerr<<"VectorL result: "
	   <<"  sumEscm = "<< vecCore::ReduceAdd(glSumEscm)
	   <<"\t  sumEkin = "<< vecCore::ReduceAdd(glSumEkin)
	   <<"\t  sumP2 = "<< vecCore::ReduceAdd(glSumP2)
	   <<"\t CPUtime = "<< glElapsed / nEvents / nReps
	   <<"\n\n\n";
  */

  //.. cleanup
  delete bulletHandler;
  delete targetHandler;
}
