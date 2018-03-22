#ifndef randomize_h
#define randomize_h 1

#include <CLHEP/Random/Randomize.h>

#define G4Random CLHEP::HepRandom

#define G4UniformRand() CLHEP::HepRandom::getTheEngine()->flat()

#endif // randomize_h 
