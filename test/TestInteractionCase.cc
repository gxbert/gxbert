//
// File:    TestInteractionCase.cpp
// Purpose: Unit tests for the vectorized class GXInteractionCase
//
// 20180611 Guilherme Lima - created

//-- ensure asserts are compiled in
#undef NDEBUG

#include "VecCore/VecCore"
#include "GXInteractionCase.hh"
#include "G4InuclParticleNames.hh"
#include "G4NucleiModel.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"

int main()
{
  using namespace gxbert;
  using namespace G4InuclParticleNames;

  // proton kinematics
  double kinEnergy = 1500.;
  double pMass = 938.272013;
  double etot = pMass + kinEnergy;
  double pmom = sqrt(etot*etot - pMass*pMass);
  LorentzVector<double> p4vec(0., 0., pmom, etot);
  GXInuclElementaryParticle<double> aProton(p4vec, proton);

  // neutron at rest
  double nMass = 939.56563;
  LorentzVector<double> n4vec(0., 0., 1., nMass);
  GXInuclElementaryParticle<double> aNeutron(n4vec, neutron);

  // build interaction case
  GXInteractionCase<double> case1(&aProton, &aNeutron);
  //std::cerr<<"IntCase(proton,neutron): "<< case1 <<"\n";

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


  //=== try a vector of GXInuclElemParticles
  using Real_v = Vc::Vector<double>;

  int vsize = vecCore::VectorSize<Real_v>();
  GXTrackHandler *handler = new GXTrackHandler(vsize);
  double pmin = sqrt((pMass + 0.5*kinEnergy)*(pMass + 0.5*kinEnergy) - pMass*pMass);
  double pmax = sqrt((pMass + 1.5*kinEnergy)*(pMass + 1.5*kinEnergy) - pMass*pMass);
  handler->GenerateRandomTracks(vsize, pmin, pmax);
  GXTrack_v track_soa = handler->GetSoATracks();
  LorentzVector<Real_v> vlorvec;
  track_soa.getFourMomentum(0,vlorvec);

  // Real_v vkinEnergy(700, 1000, 1200, 1500);
  // Real_v vpMass(pMass);
  // Real_v vetot = vpMass + vkinEnergy;
  // Real_v vpmom = Sqrt(vetot*vetot - pMass*pMass);
  // LorentzVector<Real_v> vlorvec(0., 0., vpmom, vetot);

  //std::cerr << " LorentzVector: "<< vlorvec <<"\n";

  // test default constructor -- all units based on MeV
  Index_v<Real_v> mytypes;
  Set(mytypes, 0, proton);
  Set(mytypes, 1, neutron);
  Set(mytypes, 2, photon);
  Set(mytypes, 3, deuteron);
  //Set(mytypes, 3, -29); // muonPlus -- muons are not implemented

  GXInuclParticle<Real_v> *bullets = new GXInuclElementaryParticle<Real_v>(vlorvec, mytypes);
  //std::cerr<<"\n=== GXInuclEP-bullets: "<< *dynamic_cast<GXInuclElementaryParticle<Real_v> const*>(bullets) <<"\n";

  // check re-writing
  Set(mytypes, 0, pip);
  Set(mytypes, 1, pim);
  Set(mytypes, 2, diproton);
  Set(mytypes, 3, dineutron);
  //Set(mytypes, 3, ap); // ap does not work

  GXInuclParticle<Real_v> *targets = new GXInuclElementaryParticle<Real_v>(vlorvec, mytypes);
  //std::cerr<<"\n=== GXInuclEP-targets: "<< *dynamic_cast<GXInuclElementaryParticle<Real_v> const*>(targets) <<"\n";

  GXInteractionCase<Real_v> vecCase(bullets, targets);
  //std::cerr<<"]\n*** vecCase: "<< vecCase <<"\n";

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

  //=== display result
  std::cerr<<"\n>>> GXInteractionCase tests passed.\n";
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
