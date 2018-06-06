//
// File:    TestBox.cpp
// Purpose: Unit tests for the box
//

//-- ensure asserts are compiled in
#undef NDEBUG

#include "VecCore/VecCore"
#include "GXInuclElementaryParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4NucleiModel.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"

int main()
{
  using namespace gxbert;
  using namespace G4InuclParticleNames;

  // setup kinematics
  double kinEnergy = 1500.;
  double protonMass = 938.272013;
  double etot = protonMass + kinEnergy;
  double pmom = sqrt(etot*etot - protonMass*protonMass);
  LorentzVector<double> lorvec(0., 0., pmom, etot);

  // test default constructor (all units based on MeV)
  GXInuclElementaryParticle<double> aParticle;
  aParticle.setKineticEnergy(kinEnergy);
  aParticle.setMomentumDirection( lorvec.Vect().Unit() );
  aParticle.setType(proton);

  // full constructor
  GXInuclElementaryParticle<double> aProton(lorvec, proton);

  // copy constructor
  GXInuclElementaryParticle<double> aCopy(aParticle);

  // assignment operator
  GXInuclElementaryParticle<double> assigned = aProton;

  // equality operator
  assert( aParticle == aProton);
  assert( assigned == aCopy );

  //=== test basic functionality
  assert( aProton.nucleon() );
  assert(!aProton.isMuon() );
  assert(!aProton.isPhoton());
  assert(!aProton.isElectron());
  assert(!aProton.isNeutrino());
  assert(!aProton.pion());
  assert(!aProton.antinucleon());

  //=== repeat for a neutron
  GXInuclElementaryParticle<double> aNeutron(lorvec, neutron);
  assert( aNeutron.nucleon() );
  assert(!aNeutron.isMuon() );
  assert(!aNeutron.isPhoton());
  assert(!aNeutron.isElectron());
  assert(!aNeutron.isNeutrino());
  assert(!aNeutron.pion());
  assert(!aNeutron.antinucleon());

  //=== repeat for a photon
  GXInuclElementaryParticle<double> aPhoton(lorvec, photon);
  assert(!aPhoton.nucleon() );
  assert(!aPhoton.isMuon() );
  assert( aPhoton.isPhoton());
  assert(!aPhoton.isElectron());
  assert(!aPhoton.isNeutrino());
  assert(!aPhoton.pion());
  assert(!aPhoton.antinucleon());

  //=== repeat for an electron
  GXInuclElementaryParticle<double> anElectron(lorvec, electron);
  assert(!anElectron.nucleon() );
  assert(!anElectron.isMuon() );
  assert(!anElectron.isPhoton());
  assert( anElectron.isElectron());
  assert(!anElectron.isNeutrino());
  assert(!anElectron.pion());
  assert(!anElectron.antinucleon());

  //=== repeat for a muon
  GXInuclElementaryParticle<double> aMuon(lorvec, mup);
  assert(!aMuon.nucleon() );
  assert( aMuon.isMuon() );
  assert(!aMuon.isPhoton());
  assert(!aMuon.isElectron());
  assert(!aMuon.isNeutrino());
  assert(!aMuon.pion());
  assert(!aMuon.antinucleon());

  //=== repeat for an neutrino
  GXInuclElementaryParticle<double> aNeutrino(lorvec, mnu);
  assert(!aNeutrino.nucleon() );
  assert(!aNeutrino.isMuon() );
  assert(!aNeutrino.isPhoton());
  assert(!aNeutrino.isElectron());
  assert( aNeutrino.isNeutrino());
  assert(!aNeutrino.pion());
  assert(!aNeutrino.antinucleon());

  //=== repeat for a pion
  GXInuclElementaryParticle<double> aPion(lorvec, pip);
  assert(!aPion.nucleon() );
  assert(!aPion.isMuon() );
  assert(!aPion.isPhoton());
  assert(!aPion.isElectron());
  assert(!aPion.isNeutrino());
  assert( aPion.pion());
  assert(!aPion.antinucleon());

  //=== repeat for an anti-proton
  GXInuclElementaryParticle<double> aPbar(lorvec, ap);
  assert(!aPbar.nucleon() );
  assert(!aPbar.isMuon() );
  assert(!aPbar.isPhoton());
  assert(!aPbar.isElectron());
  assert(!aPbar.isNeutrino());
  assert(!aPbar.pion());
  assert( aPbar.antinucleon());

  //=== try a vector of GXInuclElemParticles
  using Real_v = Vc::Vector<double>;

  int vsize = vecCore::VectorSize<Real_v>();
  GXTrackHandler *handler = new GXTrackHandler(vsize);
  double pmin = sqrt((protonMass + 0.5*kinEnergy)*(protonMass + 0.5*kinEnergy) - protonMass*protonMass);
  double pmax = sqrt((protonMass + 1.5*kinEnergy)*(protonMass + 1.5*kinEnergy) - protonMass*protonMass);
  handler->GenerateRandomTracks(vsize, pmin, pmax);
  GXTrack_v track_soa = handler->GetSoATracks();
  LorentzVector<Real_v> vlorvec;
  track_soa.getFourMomentum(0,vlorvec);

  // Real_v vkinEnergy(700, 1000, 1200, 1500);
  // Real_v vprotonMass(protonMass);
  // Real_v vetot = vprotonMass + vkinEnergy;
  // Real_v vpmom = Sqrt(vetot*vetot - protonMass*protonMass);
  // LorentzVector<Real_v> vlorvec(0., 0., vpmom, vetot);

  //std::cerr << " LorentzVector: "<< vlorvec <<"\n";

  // test default constructor -- all units based on MeV
  Index_v<Real_v> mytypes;
  Set(mytypes, 0, proton);
  Set(mytypes, 1, neutron);
  Set(mytypes, 2, photon);
  Set(mytypes, 3, deuteron);
  //Set(mytypes, 3, -29); // muonPlus -- muons are not implemented

  GXInuclElementaryParticle<Real_v> tracks(vlorvec, mytypes);
  std::cerr<<"=== GXInuclElemParticles: "<< tracks <<"\n";
  std::cerr<< vsize <<" tracks: "<< tracks.type()
	   <<"\n kinE="<< tracks.getKineticEnergy()
	   <<"\n totE="<< tracks.getTotalEnergy()
	   <<"\n nucleon:"<< tracks.nucleon()
	   <<"\n pion:"<< tracks.pion()
	   <<"\n photon:"<< tracks.isPhoton()
	   <<"\n baryon:"<< tracks.baryon()
	   <<"\n strange:"<< tracks.getStrangeness()
	   <<"\n quasi_deutron(): "<< tracks.quasi_deutron()
	   <<"\n";

  // check re-writint
  Set(mytypes, 0, pip);
  Set(mytypes, 1, pim);
  Set(mytypes, 2, diproton);
  Set(mytypes, 3, dineutron);
  //Set(mytypes, 3, ap); // ap does not work

  tracks.fill(vlorvec, mytypes);
  std::cerr<<"=== GXInuclElemParticles: "<< tracks <<"\n";
  std::cerr<< vsize <<" tracks: "<< tracks.type()
	   <<"\n kinE="<< tracks.getKineticEnergy()
	   <<"\n totE="<< tracks.getTotalEnergy()
	   <<"\n nucleon:"<< tracks.nucleon()
	   <<"\n pion:"<< tracks.pion()
	   <<"\n photon:"<< tracks.isPhoton()
	   <<"\n baryon:"<< tracks.baryon()
	   <<"\n strange:"<< tracks.getStrangeness()
	   <<"\n quasi_deutron(): "<< tracks.quasi_deutron()
	   <<"\n";

  //=== display result
  std::cerr<<">>> GXInuclElementaryParticle tests passed.\n";
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
