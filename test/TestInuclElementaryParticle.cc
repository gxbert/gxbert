//
// File:    TestBox.cpp
// Purpose: Unit tests for the box
//

//-- ensure asserts are compiled in
#undef NDEBUG

#include "GXInuclElementaryParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4NucleiModel.hh"

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
}
