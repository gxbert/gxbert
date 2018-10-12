//
// File: ElemParticleColliderBenchmark.cc
//
// Purpose: Benchmark for G4/GX ElementaryParticleCollider classes
//
// 2018-05-25 Guilherme lima  Adapted from InuclColliderBenchmark.cc

#include "VecCore/VecCore"
#include "timer.h"
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
#include "Randomize.hh"

//using namespace ::vecCore;
using namespace gxbert;
using std::cerr;

//.. default globals
constexpr double pmass = 0.938272013;  // proton mass in GeV
size_t nReps     = 1;
size_t nEvents   = 1024; // 1024*64;
double kinEnergy = 1.500; // in GeV

int debugLevel   = 0;

bool isOutputInvalid(G4InuclParticle const &bullet, G4CollisionOutput const &output)
{
  // Quantities necessary for conditional block below
  G4int npart = output.numberOfOutgoingParticles();
  G4int nfrag = output.numberOfOutgoingNuclei();

  const GXParticleDefinition* out0 = (npart == 0) ? 0 : output.getOutgoingParticles()[0].getDefinition();
  const GXParticleDefinition* out1 = (npart == 0) ? 0 : output.getOutgoingParticles()[1].getDefinition();

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

  // return ( ((npart != 0) &&
    // 	    //#ifdef G4CASCADE_COULOMB_DEV
    // 	    //(coulombBarrierViolation() && npart+nfrag > 2)
    // 	    //#else
    // 	      (npart+nfrag < 3 && (out0 == bullet.getDefinition() || out1 == bullet.getDefinition()))
    // 	    //#endif
    //          )
    // 	     //#ifndef G4CASCADE_SKIP_ECONS
    // 	     //|| (!balance->okay())
    // 	     //#endif
    // 	 );

  bool result = (npart==0);

  if(debugLevel>1) {
    std::cout<<"isOutputInvalid/G4: npart="<< npart <<", nfrag="<< nfrag
	     <<", out0="<< out0->GetParticleName()
	     <<", out1="<< out1->GetParticleName()
	     <<", bullet="<< bullet.getDefinition()->GetParticleName()
	     <<" --> return "<< result <<"\n";
  }

  return result;
}

template <typename T>
bool isOutputInvalid(GXInuclParticle<T> const& bullet, GXCollisionOutput<T> const& output)
{
  // Quantities necessary for conditional block below
  int npart = output.numberOfOutgoingParticles();
  int nfrag = output.numberOfOutgoingNuclei();

  Index_v<T> types = output.getOutgoingParticles().begin()->type();
  //const GXParticleDefinition* firstOut = (npart == 0) ? 0 : output.getOutgoingParticles().begin()->getDefinition();

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

  Index_v<T> btypes = dynamic_cast<GXInuclElementaryParticle<T> const*>(&bullet)->type();
  // GXInuclElementaryParticle<T> const* epBullet =
  // GXParticleDefinition const* bDef0 = (epBullet ? getDefinition(Get(epBullet->type(),0)) : 0);

  // return ( ((npart != 0) &&
  // 	    //#ifdef G4CASCADE_COULOMB_DEV
  // 	    //(coulombBarrierViolation() && npart+nfrag > 2)
  // 	    //#else
  // 	    (npart+nfrag < 3 && MaskFull(types == btypes))
  // 	    //#endif
  //            )
  // 	     //#ifndef G4CASCADE_SKIP_ECONS
  // 	     //|| (!balance->okay())
  // 	     //#endif
  // 	 );

  bool result = (npart==0);
  if(debugLevel>1) {
    std::cout<<"isOutputInvalid/G4: npart="<< npart <<", nfrag="<< nfrag <<", firstTypes="<< types <<", bullet.types="<< btypes
	     <<" --> return "<< result <<"\n";
  }

  return result;
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
    for(size_t i = 0; i < nEvents; ++i) {
      soaBullets.getFourMomentum(i, lorvec);
      bullet.fill(lorvec, G4InuclParticleNames::proton);

      soaTargets.getFourMomentum(i, lorvec);
      target.fill(lorvec, G4InuclParticleNames::neutron);

      if (debugLevel>1) {
	std::cerr<<"\n=== G4InuclEP-bullet: "<< bullet;
	std::cerr<<"\n=== G4InuclEP-target: "<< target <<"\n\n";
      }

      numberOfTries = 0;
      do {   			// create inelastic interaction
	output.reset();
	collider.collide(&bullet, &target, output);
	numberOfTries++;
      } while ( numberOfTries <= 20 && isOutputInvalid(bullet, output) );

      // update counters
      if (debugLevel > 0) {
	cerr<<"***** Event "<< i <<" Ntries="<< numberOfTries <<" Final state: "<< output.numberOfOutgoingParticles()
	    <<" hadrons + "<< output.numberOfOutgoingNuclei() <<" nuclei\n";
      }

      nHadrons += output.numberOfOutgoingParticles();
      const std::vector<G4InuclElementaryParticle>& outParticles = output.getOutgoingParticles();
      std::vector<G4InuclElementaryParticle>::const_iterator ipart = outParticles.begin(), iend = outParticles.end();
      for( ; ipart != iend; ++ipart) {
       	if (ipart->type() == G4InuclParticleNames::neutron) ++nNeutrons;
       	else if (ipart->type() == G4InuclParticleNames::proton) ++nProtons;
       	else if (ipart->type() == G4InuclParticleNames::photon) ++nPhotons;
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
void RunGXElementaryParticleCollider(const char* testname, GXTrack_v const& soaBullets, GXTrack_v const& soaTargets,
				     GXCollisionOutput<Real_v> &output)
{
  int nNeutrons = 0;
  int nProtons = 0;
  int nPhotons = 0;
  int nOthers = 0;

  // build a vectorized EP collider
  int vsize = vecCore::VectorSize<Real_v>();
  GXElementaryParticleCollider<Real_v> collider;
  collider.setVerboseLevel(debugLevel);

  GXInuclElementaryParticle<Real_v> bullets, targets;
  LorentzVector<Real_v> lorvec;
  Timer<nanoseconds> timer;
  //vecCore::Index_v<Real_v> counter(0);
  int numberOfTries = 0;

  // allocate buffer for multiplicity array
  // size_t intsize(vecCore::Index_v<Real_v>);
  // vecCore::Index_v<Real_v>* pMult = (vecCore::Index_v<Real_v>*)_mm_malloc(nEvents*intsize, 64);
  //long pMult[nEvents];

  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t iv = 0, i = 0; iv < nEvents; ++i, iv += vsize) {
      soaBullets.getFourMomentum(iv, lorvec);
      bullets.fill(lorvec, G4InuclParticleNames::proton);

      soaTargets.getFourMomentum(iv, lorvec);
      targets.fill(lorvec, G4InuclParticleNames::neutron);

      GXInuclElementaryParticle<Real_v> const* pbullets = dynamic_cast<GXInuclElementaryParticle<Real_v> const*>(&bullets);
      GXInuclElementaryParticle<Real_v> const* ptargets = dynamic_cast<GXInuclElementaryParticle<Real_v> const*>(&targets);
      if (debugLevel>1) {
	std::cerr<<"\n=== GXInuclEP-bullets: "<< *pbullets;
	std::cerr<<"\n=== GXInuclEP-targets: "<< *ptargets <<"\n\n";
      }

      /*
      //std::cerr<<"useEPCollider<Real_v>(bullets,targtes): "<< collider.useEPCollider(pbullets, ptargets) <<"\n";
      // assert(collider.useEPCollider(pbullets, ptargets));
      // assert(collider.useEPCollider(ptargets, pbullets));

      auto initStateVec = pbullets->type() * ptargets->type();
      Real_v ekinVec = bullets.getKineticEnergy();

      auto multipl = collider.generateMultiplicity(initStateVec, ekinVec);
      counter += multipl;
      pMult[i] = multipl;

      collider.generateOutgoingPartTypes(initStateVec, multipl, ekinVec);
      collider.fillOutgoingMasses();
      */

      // only needed for vectorized mode
      //collider.rebasketizeByMultiplicity(nEvents, pMult, &bullets, &targets);

      if (debugLevel > 1) cerr<<"=== Generating cascade attemp "<< numberOfTries <<"\n";

      numberOfTries = 0;
      do {   			// we try to create inelastic interaction
	output.reset();
	collider.collide(pbullets, ptargets, output);
	numberOfTries++;
      } while ( numberOfTries <= 20 && isOutputInvalid<Real_v>(*pbullets, output) );

      // update counters
      if (debugLevel>0) {
	cerr<<"***** Event "<< i <<" Ntries="<< numberOfTries
	    <<" Final state: "<< output.numberOfOutgoingParticles()
	    <<" hadrons + 0 nuclei\n"; //<< output.numberOfOutgoingNuclei() <<" nuclei\n";
      }

      /*
      auto initStateVec = soaBullets.type() * soaTargets.type();
      Real_v ekinVec = soaBullets.getKineticEnergy();
      nHadrons += collider.generateMultiplicity(initStateVec, ekinVec);
      std::cerr<<" nHadrons: "<< nHadrons <<"\n";
      */

      const std::vector<GXInuclElementaryParticle<Real_v>>& outParticles = output.getOutgoingParticles();
      typename std::vector<const GXInuclElementaryParticle<Real_v>>::const_iterator ipart = outParticles.begin();
      typename std::vector<const GXInuclElementaryParticle<Real_v>>::const_iterator iend  = outParticles.end();
      for( ; ipart != iend; ++ipart) {
	Index_v<Real_v> types = ipart->type();
	for(size_t j = 0; j < VectorSize<Real_v>(); ++j) {
	  int jtype = Get(types, j);
	  if (jtype == G4InuclParticleNames::neutron) ++nNeutrons; //??? suspended...
	  else if (jtype == G4InuclParticleNames::proton) ++nProtons; //??? suspended...
	  else if (jtype == G4InuclParticleNames::photon) ++nPhotons; //??? suspended...
	  else ++nOthers;
	}
      }
    }

    //_mm_free(pMult);

  }
  double elapsed = timer.Elapsed();

  std::cerr<< testname <<" results: "
	   <<"  nHadrons="<< (nNeutrons + nProtons + nOthers)
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
  bulletHandler->SetMass( 0.93827203 );  // proton mass in GeV
  double posdir[6] = {0., 0., 0., 0., 0., 1.};
  bulletHandler->GenerateTracksAlongSameDirection(nEvents, posdir, pmom, pmom);
  GXTrack_v &soaBullets = bulletHandler->GetSoATracks();
  printf("Number of tracks = %d at Ekin=%f GeV and pMom=%f GeV\n", soaBullets.size, kinEnergy, pmom);

  //.. prepare targets at ~rest
  GXTrackHandler *targetHandler = new GXTrackHandler(nEvents);
  targetHandler->SetMass( 0.93956536 );  // neutron mass in GeV
  targetHandler->GenerateTracksAlongSameDirection(nEvents, posdir, 1.0e-6, 1.0e-6);
  GXTrack_v &soaTargets = targetHandler->GetSoATracks();
  // for(int i=0; i<nEvents; ++i) {
  //   soaTargets.E[i] = 0.; // kinEnergy
  // }

  // Using function calls for benchmarks
  //CLHEP::HepRandomEngine* theEngine = G4Random::getTheEngine();
  CLHEP::HepRandom::setTheSeed(17263543);
  G4CollisionOutput output;
  RunG4ElementaryParticleCollider(soaBullets, soaTargets, output);

  // benchmarks templated on types
  //theEngine->setSeeds(seeds);
  CLHEP::HepRandom::setTheSeed(17263543);
  GXCollisionOutput<double> scalarOutput;
  RunGXElementaryParticleCollider<double>("double", soaBullets, soaTargets, scalarOutput);

  //theEngine->setSeeds(seeds);
  CLHEP::HepRandom::setTheSeed(17263543);
  GXCollisionOutput<Real_v> vectorOutput;
  RunGXElementaryParticleCollider<Real_v>("Real_v", soaBullets, soaTargets, vectorOutput);

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
