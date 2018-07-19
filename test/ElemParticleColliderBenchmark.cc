//
// File: ElemParticleColliderBenchmark.cc
//
// Purpose: Benchmark for G4/GX ElementaryParticleCollider classes
//
// 2018-05-25 Guilherme lima  Adapted from InuclColliderBenchmark.cc

#include "VecCore/VecCore"
#include "VecCore/Timer.h"
#include "GXTrack.hh"
#include "GXTrackHandler.hh"
//#include "GXInuclCollider.hh"

#include "G4InuclElementaryParticle.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4CollisionOutput.hh"

#include "LorentzVector.hh"
#include "GXElementaryParticleCollider.hh"
#include "GXParticleDefinition.hh"
#include "GXCollisionOutput.hh"

#include "GXProton.hh"
#include "GXNeutron.hh"
#include "G4InuclParticleNames.hh"

#include <iostream>

//using namespace ::vecCore;
using namespace gxbert;
using std::cerr;

//.. default globals
constexpr double pmass = 0.938272013;  // proton mass in GeV
size_t nReps     = 1;
size_t nEvents   = 16; // 1024*64;
double kinEnergy = 1.500; // in GeV

int debugLevel   = 0;

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


//
//================================= Reference: G4InuclCollider functionalities ==================
//
void RunG4ElementaryParticleCollider(GXTrack_v const& soaBullets, GXTrack_v const& soaTargets, G4CollisionOutput &output)
{
  int nHadrons = 0;
  int nNeutrons = 0;
  int nProtons = 0;
  int nPhotons = 0;
  int nOthers = 0;

  G4ElementaryParticleCollider collider;
  collider.setVerboseLevel(debugLevel);

  G4InuclElementaryParticle bullet, target;
  G4LorentzVector lorvec;
  Timer<nanoseconds> timer;
  timer.Start();
  int numberOfTries = 0;
  for(size_t irep = 0; irep < nReps; ++irep) {
    nHadrons = nNeutrons = nProtons = nPhotons = nOthers = 0;
    for(size_t i = 0; i < nEvents; ++i) {
      soaBullets.getFourMomentum(i, lorvec);
      bullet.fill(lorvec, G4InuclParticleNames::proton);

      soaTargets.getFourMomentum(i, lorvec);
      target.fill(lorvec, G4InuclParticleNames::proton);

      numberOfTries = 0;
      do {   			// we try to create inelastic interaction
	if (debugLevel > 1) cerr<<"=== Generating cascade attempt " << numberOfTries <<"\n";

	output.reset();
	collider.collide(&bullet, &target, output);

	numberOfTries++;
      } while ( numberOfTries <= 20 && isOutputInvalid(bullet, output) );

      // update counters
      cerr<<"***** Event "<< i <<" Ntries="<< numberOfTries <<" Final state: "<< output.numberOfOutgoingParticles() <<" hadrons + "<< output.numberOfOutgoingNuclei() <<" nuclei\n";
      nHadrons += output.numberOfOutgoingParticles();
      const std::vector<G4InuclElementaryParticle>& outParticles = output.getOutgoingParticles();
      std::vector<G4InuclElementaryParticle>::const_iterator ipart = outParticles.begin(), iend = outParticles.end();
      for( ; ipart != iend; ++ipart) {
       	if (ipart->type() == G4InuclParticleNames::neutron) ++nNeutrons; //??? suspended...
       	else if (ipart->type() == G4InuclParticleNames::proton) ++nProtons; //??? suspended...
       	else if (ipart->type() == G4InuclParticleNames::photon) ++nPhotons; //??? suspended...
	else ++nOthers;
      }
    }
  }
  double g4Elapsed = timer.Elapsed();

  std::cerr<<"GXBert results: "
	   <<"  nHadrons="<< nHadrons
	   <<"\tnNeutrons="<< nNeutrons
	   <<"\tnProtons="<< nProtons
	   <<"\tnPhotons="<< nPhotons
	   <<"\tnOthers="<< nOthers
	   <<"\tCPUtime="<< g4Elapsed / nEvents / nReps
	   <<"\n";
}


//
//=====================================  VECTORIZED testing version  ========================
//
template <typename Real_v>
void RunGXElementaryParticleCollider(const char* testname, GXTrack_v const& soaBullets, GXTrack_v const& soaTargets, GXCollisionOutput &output)
{
  int nHadrons = 0;
  int nNeutrons = 0;
  int nProtons = 0;
  int nPhotons = 0;
  int nOthers = 0;

  int vsize = vecCore::VectorSize<Real_v>();
  GXElementaryParticleCollider<Real_v> collider;
  collider.setVerboseLevel(debugLevel);

  GXInuclElementaryParticle<Real_v> bullet, target;
  LorentzVector<Real_v> lorvec;
  Timer<nanoseconds> timer;
  timer.Start();
  int numberOfTries = 0;
  for(size_t irep = 0; irep < nReps; ++irep) {
    nHadrons = nNeutrons = nProtons = nPhotons = nOthers = 0;
    for(size_t i = 0; i < nEvents; i += vsize) {
      soaTargets.getFourMomentum(i, lorvec);
      bullet.fill(lorvec, Real_v(G4InuclParticleNames::proton));

      soaBullets.getFourMomentum(i, lorvec);
      target.fill(lorvec, Real_v(G4InuclParticleNames::proton));

      numberOfTries = 0;

      if (debugLevel > 1) cerr<<"=== Generating cascade attemp "<< numberOfTries <<"\n";
      output.reset();
      collider.collide(&bullet, &target, output);
      ++numberOfTries;

      // update counters
      cerr<<"***** Event "<< i <<" Ntries="<< numberOfTries <<"\n";
      //<<" Final state: "<< output.numberOfOutgoingParticles()
      //<<" hadrons + "<< output.numberOfOutgoingNuclei() <<" nuclei\n";
      //nHadrons += output.numberOfOutgoingParticles();
      //const std::vector<G4InuclElementaryParticle>& outParticles = output.getOutgoingParticles();
      //std::vector<const G4InuclElementaryParticle>::const_iterator ipart = outParticles.begin(), iend = outParticles.end();
      //for( ; ipart != iend; ++ipart) {
      //  if (ipart->type() == G4InuclParticleNames::neutron) ++nNeutrons; //??? suspended...
      //  else if (ipart->type() == G4InuclParticleNames::proton) ++nProtons; //??? suspended...
      //  else if (ipart->type() == G4InuclParticleNames::photon) ++nPhotons; //??? suspended...
      //  else ++nOthers;
      // }
    }
  }
  double elapsed = timer.Elapsed();

  std::cerr<< testname <<" results: "
	   <<"  nHadrons="<< nHadrons
	   <<"\tnNeutrons="<< nNeutrons
	   <<"\tnProtons="<< nProtons
	   <<"\tnPhotons="<< nPhotons
	   <<"\tnOthers="<< nOthers
	   <<"\tCPUtime="<< elapsed / nEvents / nReps
	   <<"\n";
}


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

  // Note 1: both MeV and GeV can be used, but it needs to be consistent!
  // Note 2: Bertini algorithms expect units in GeV

  //.. prepare bullets
  GXTrackHandler *bulletHandler = new GXTrackHandler(nEvents); 
  bulletHandler->SetMass( 0.938272013 );  // proton mass in GeV
  double posdir[6] = {0., 0., 0., 0., 0., 1.};
  bulletHandler->GenerateTracksAlongSameDirection(nEvents, posdir, pmom, pmom);
  GXTrack_v &soaBullets = bulletHandler->GetSoATracks();
  printf("Number of tracks = %d at Ekin=%f GeV and pMom=%f GeV\n", soaBullets.size, kinEnergy, pmom);

  //.. prepare targets at ~rest
  GXTrackHandler *targetHandler = new GXTrackHandler(nEvents);
  targetHandler->SetMass( 0.938272013 );  // proton mass in GeV
  targetHandler->GenerateTracksAlongSameDirection(nEvents, posdir, 0.001, 0.001);
  GXTrack_v &soaTargets = targetHandler->GetSoATracks();
  // for(int i=0; i<nEvents; ++i) {
  //   soaTargets.E[i] = 0.; // kinEnergy
  // }

  // Using function calls for benchmarks
  G4CollisionOutput output;
  RunG4ElementaryParticleCollider(soaBullets, soaTargets, output);

  // benchmarks templated on types
  GXCollisionOutput gxoutput;
  RunGXElementaryParticleCollider<double>("double", soaBullets, soaTargets, gxoutput);
  //RunGXElementaryParticleCollider<vecCore::backend::VcVector::Double_v>("Double_v", soaBullets, soaTargets);

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
  return 0;
}
