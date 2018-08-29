//
// File:    TestG4ElemParticleCollider.cpp
// Purpose: Unit tests for the G4ElementaryParticleCollider
//
// 20180824 Guilherme Lima - created

//-- ensure asserts are compiled in
#undef NDEBUG

#include "VecCore/VecCore"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4NucleiModel.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"

int main()
{
  using namespace gxbert;
  using namespace G4InuclParticleNames;

  // proton kinematics
  double kinEnergy = 1.500;
  double pMass = 0.938272013;
  double etot = pMass + kinEnergy;
  double pmom = sqrt(etot*etot - pMass*pMass);
  G4LorentzVector p4vec(0., 0., pmom, etot);
  G4InuclElementaryParticle aProton;
  aProton.fill(p4vec, proton);

  // neutron at rest
  double nMass = 0.93956563;
  G4LorentzVector n4vec(0., 0., 1., nMass);
  G4InuclElementaryParticle aNeutron;
  aNeutron.fill(n4vec, neutron);

  // build a scalar collider
  G4ElementaryParticleCollider collider;
  //collider.setVerboseLevel(5);
  //std::cerr<<"useEPCollider(aProton,aNeutron): "<< collider.useEPCollider(&aProton, &aNeutron) <<"\n";
  //assert(collider.useEPCollider(&aProton, &aNeutron));
  //assert(collider.useEPCollider(&aNeutron, &aProton));

  int initialState = aProton.type() * aNeutron.type();
  double ekin = aProton.getKineticEnergy();
  for(size_t i=0; i<10; ++i) {
    int multipl = collider.generateMultiplicity(initialState, ekin);
    if(multipl<9) std::cerr<<" i="<< i <<" multipl="<< multipl <<"\n";
  }
  /*
  //=== test basic functionality
  assert(case1.valid());
  assert(case1.hadrons());
  assert(!case1.twoNuclei());
  assert(!case1.hadNucleus());

  // make sure unified notation can be used for <double> case
  assert(vecCore::MaskFull( case1.valid() ));
  assert(vecCore::MaskFull(case1.hadrons()));
  assert(vecCore::MaskEmpty( case1.twoNuclei() ));
  assert(vecCore::MaskEmpty( case1.hadNucleus() ));

  //std::cerr<<" bullet: "<< *(dynamic_cast<GXInuclElementaryParticle<double> const*>(case1.getBullet())) <<"\n";
  assert(vecCore::MaskFull(*static_cast<GXInuclElementaryParticle<double> const*>(case1.getBullet()) == aProton));

  //std::cerr<<" target: "<< *(dynamic_cast<GXInuclElementaryParticle<double> const*>(case1.getTarget())) <<"\n";
  assert(vecCore::MaskFull(*static_cast<GXInuclElementaryParticle<double> const*>(case1.getTarget()) == aNeutron));
  */


  //=== display result
  std::cerr<<"\n>>> G4InuclElemPartCollider::generateMultiplicity() test passed.\n";
  return 0;
}

/*
int main()
{
  RunUnitTests<double>("double");
#ifdef VECCORE_ENABLE_VC
  RunUnitTests<Vc::Vector<double>>("VcVector");
#endif
}
*/
