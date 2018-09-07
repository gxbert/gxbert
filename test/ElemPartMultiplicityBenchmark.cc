//
// File: ElemPartMultiplicityBenchmark.cc
//
// Purpose: Benchmark for G4/GX ElementaryParticleCollider function generateMultiplicity()
//
// Note:    Requires this method to be made public in G*ElementaryParticleCollider class
//
// 2018-08-22 Guilherme lima  Adapted from ElemParticleColliderBenchmark.cc

// this is needed in order to make generateMultiplicity public
#undef NDEBUG
#include "VecCore/VecCore"
#include "VecCore/Timer.h"
#include "GXTrack.hh"
#include "GXTrackHandler.hh"

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
size_t nReps     = 100;
size_t nEvents   = (1 << 20); // 1024*64;
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
void RunG4ElemParticleMultiplicity(GXTrack_v const& soaBullets, GXTrack_v const& soaTargets, G4CollisionOutput &output)
{
  int counter = 0;

  G4ElementaryParticleCollider collider;
  collider.setVerboseLevel(debugLevel);

  G4InuclElementaryParticle bullet, target;
  G4LorentzVector lorvec;
  Timer<nanoseconds> timer;
  timer.Start();
  //int numberOfTries = 0;
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < nEvents; ++i) {
      soaBullets.getFourMomentum(i, lorvec);
      bullet.fill(lorvec, G4InuclParticleNames::proton);

      soaTargets.getFourMomentum(i, lorvec);
      target.fill(lorvec, G4InuclParticleNames::proton);

      int is = bullet.type() * target.type();
      double ekin = bullet.getKineticEnergy();
      int multiplicity = collider.generateMultiplicity(is, ekin);
      counter += multiplicity;
      //if (debugLevel > 3) std::cerr << " is=" << is <<" kinE="<< ekin <<", multiplicity="<< multiplicity <<"\n";
    }
  }
  double g4Elapsed = timer.Elapsed();

  std::cerr<<"GXBert results: "
	   <<"  integrated multiplicity="<< counter
	   <<"\tCPUtime="<< g4Elapsed / nEvents / nReps
	   <<"\n";
}


//
//=====================================  VECTORIZED testing version  ========================
//
template <typename Real_v>
void RunGXElemParticleMultiplicity(const char* testname, GXTrack_v const& soaBullets, GXTrack_v const& soaTargets, GXCollisionOutput &output)
{
  const size_t vsize = vecCore::VectorSize<Real_v>();
  //using Int_v = vecCore::VcSimdArray<vsize>;

  GXElementaryParticleCollider<Real_v> collider;
  collider.setVerboseLevel(debugLevel);

  GXInuclElementaryParticle<Real_v> bullet, target;
  LorentzVector<Real_v> lorvec;
  vecCore::Index_v<Real_v> counter(0);
  Timer<nanoseconds> timer;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < nEvents; i += vsize) {
      soaBullets.getFourMomentum(i, lorvec);
      bullet.fill(lorvec, G4InuclParticleNames::proton);

      soaTargets.getFourMomentum(i, lorvec);
      target.fill(lorvec, G4InuclParticleNames::neutron);

      auto initialState = bullet.type() * target.type();
      Real_v kinEnergy = bullet.getKineticEnergy();
      auto multiplicity = collider.generateMultiplicity(initialState, kinEnergy);

      // update counters
      counter += multiplicity;

      if (debugLevel > 2) {
	cerr <<"***** Track "<< i
	     <<" initState="<< initialState
	     <<", kinEnergy="<< kinEnergy
	     <<", multiplicity = "<< multiplicity <<"\n";
      }

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
	   <<"  integrated multiplicity="<< vecCore::ReduceAdd(counter)
	   <<"\tCPUtime="<< elapsed / nEvents / nReps
	   <<"\n";
}


int main(int argc, char* argv[]) {

  std::cerr<<"\n*** Usage: "<< argv[0] <<" [debugLevel [nEvents [nReps [kinEnergy [debugLevel]..]\n\n";

  if (argc >= 2) nEvents = atoi(argv[1]);
  if (argc >= 3) nReps = atoi(argv[2]);
  if (argc >= 4) kinEnergy  = atof(argv[3]); // energy in GeV
  if (argc >= 5) debugLevel = atoi(argv[4]);

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
  bulletHandler->GenerateTracksAlongSameDirection(nEvents, posdir, 0.9*pmom, 1.1*pmom);
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
  RunG4ElemParticleMultiplicity(soaBullets, soaTargets, output);

  // benchmarks templated on types
  GXCollisionOutput gxoutput;
  RunGXElemParticleMultiplicity<double>("double", soaBullets, soaTargets, gxoutput);
  RunGXElemParticleMultiplicity<Real_v>("Real_v", soaBullets, soaTargets, gxoutput);

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
