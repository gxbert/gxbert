//
// File:    TestInteractionCase.cpp
// Purpose: Unit tests for the vectorized class GXInteractionCase
//
// 20180611 Guilherme Lima - created

//-- ensure asserts are compiled in
#undef NDEBUG

#include "VecCore/VecCore"
#include "GXElementaryParticleCollider.hh"
#include "G4InuclParticleNames.hh"
#include "G4NucleiModel.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"

constexpr double MeV = 0.001;

int main()
{
  using namespace gxbert;
  using namespace G4InuclParticleNames;

  // proton kinematics
  double kinEnergy = 1500. * MeV;
  double pMass = 938.272013 * MeV;
  double etot = pMass + kinEnergy;
  double pmom = sqrt(etot*etot - pMass*pMass);
  LorentzVector<double> p4vec(0., 0., pmom, etot);
  GXInuclElementaryParticle<double> aProton(p4vec, proton);

  // neutron at rest
  double nMass = 939.56563 * MeV;
  LorentzVector<double> n4vec(0., 0., 1., nMass);
  GXInuclElementaryParticle<double> aNeutron(n4vec, neutron);

  // build a scalar collider
  GXElementaryParticleCollider<double> collider;
  collider.setVerboseLevel(5);
  std::cerr<<"useEPCollider(aProton,aNeutron): "<< collider.useEPCollider(&aProton, &aNeutron) <<"\n";
  assert(collider.useEPCollider(&aProton, &aNeutron));
  assert(collider.useEPCollider(&aNeutron, &aProton));

  int initialState = aProton.type() * aNeutron.type();
  double ekin = aProton.getKineticEnergy();
  for(size_t i=0; i<10; ++i) {
    auto multipl = collider.generateMultiplicity(initialState, ekin);
    std::cerr<<" i="<< i <<" multipl="<< multipl <<"\n";
    collider.generateOutgoingPartTypes(initialState, multipl, ekin);
    collider.fillOutgoingMasses();
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

  //=== try a vector of GXInuclElemParticles
  //using Real_v = Vc::Vector<double>;
  using Real_v = Double_v;

  int vsize = vecCore::VectorSize<Real_v>();
  GXTrackHandler *handler = new GXTrackHandler(vsize);
  double pmin = sqrt((pMass + 0.5*kinEnergy)*(pMass + 0.5*kinEnergy) - pMass*pMass);
  double pmax = sqrt((pMass + 1.5*kinEnergy)*(pMass + 1.5*kinEnergy) - pMass*pMass);
  double posdir[6] = {0., 0., 0., 0., 0., 1.};
  handler->GenerateTracksAlongSameDirection(vsize, posdir, pmin, pmax);

  GXTrack_v track_soa = handler->GetSoATracks();
  std::cerr<<"GTrack_v.size = "<< track_soa.size <<"\n";
  LorentzVector<Real_v> vlorvec = track_soa.getFourMomentum<Real_v>(0);
  std::cerr<<"\nvlorvec="<< vlorvec <<" dotProd="<< vlorvec.Mag() <<"\n";
  //track_soa.getFourMomentum(0,vlorvec);

  // Real_v vkinEnergy(700, 1000, 1200, 1500);
  // Real_v vpMass(pMass);
  // Real_v vetot = vpMass + vkinEnergy;
  // Real_v vpmom = Sqrt(vetot*vetot - pMass*pMass);
  // LorentzVector<Real_v> vlorvec(0., 0., vpmom, vetot);

  //std::cerr << " LorentzVector: "<< vlorvec <<"\n";

  // test default constructor -- all units based on GeV
  Index_v<Real_v> btypes(proton);

  GXInuclElementaryParticle<Real_v> *bullets = new GXInuclElementaryParticle<Real_v>(vlorvec, btypes);
  //std::cerr<<"\n=== GXInuclEP-bullets: "<< *dynamic_cast<GXInuclElementaryParticle<Real_v> const*>(bullets) <<"\n";

  // check re-writing
  Index_v<Real_v> ttypes(neutron);

  GXInuclElementaryParticle<Real_v> *targets = new GXInuclElementaryParticle<Real_v>(vlorvec, ttypes);
  //std::cerr<<"\n=== GXInuclEP-targets: "<< *dynamic_cast<GXInuclElementaryParticle<Real_v> const*>(targets) <<"\n";

  // build a vectorized EP collider
  GXElementaryParticleCollider<Real_v> vecCollider;
  std::cerr<<"useEPCollider<Real_v>(bullets,targtes): "<< vecCollider.useEPCollider(bullets,targets) <<"\n";
  assert(vecCollider.useEPCollider(bullets, targets));
  assert(vecCollider.useEPCollider(targets, bullets));

  auto initStateVec = bullets->type() * targets->type();
  Real_v ekinVec = bullets->getKineticEnergy();
  auto multipl = vecCollider.generateMultiplicity(initStateVec, ekinVec);
  std::cerr<<" multipl="<< multipl <<"\n";
  for(size_t i=2; i<10; ++i) {
    // make multiplicity vector uniform
    multipl = vecCore::Index_v<Real_v>(i);
    vecCollider.generateOutgoingPartTypes(initStateVec, multipl, ekinVec);
    vecCollider.fillOutgoingMasses();
  }

  /*
  // boolean test cases
  assert(vecCore::MaskFull( vecCase.valid() ));
  assert(vecCore::MaskEmpty( vecCase.twoNuclei() ));
  assert(vecCore::MaskEmpty( vecCase.hadNucleus() ));
  //std::cerr<< "hadrons():"<< vecCase.hadrons() <<"\n";
  assert(vecCore::MaskFull( vecCase.hadrons() ));

  //std::cerr<<" code:"<< vecCase.code() <<"\n";
  //assert(vecCore::MaskFull( vecCase.code() ));

  //delete bullets;  // these deletes cause "Illegal instruction: 4"  on Clang/OSX
  //delete targets;
  */

  //=== display result
  std::cerr<<"\n>>> GXElemParticleCollider tests passed.\n";
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
