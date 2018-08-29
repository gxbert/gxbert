//
// Description: class containing parametrized momentum distributions
//              in the CM for nucleon/nucleon >= 4-body final states
//
#include "GXNuclNucl4BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const double pqprC[2][4][4] = {
    { { 1.6208, -0.2009,  0.0126, -0.0002},
      {-4.3139,  1.3641, -0.0835,  0.0014},
      { 12.291,  -3.403,   0.186, -0.0024},
      {-15.288,  3.8559, -0.2004,  0.0022} },
    { {1.2419,  -0.244,  0.0157, -0.0003},
      {-4.3633,  1.3158, -0.0826,  0.0014},
      {13.743, -3.5691,  0.2143, -0.0034},
      {-18.592,  4.3867, -0.2585,  0.0039} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const double psC[2][3] = {
    { 0.6296, 0.1787, -0.0026 }, { 0.8381, 0.0086, 0.0033 }
  };
}

// Constructor passes arrays to templated base class

gxbert::GXNuclNucl4BodyMomDst::GXNuclNucl4BodyMomDst(int verbose)
  : GXVMultiBodyMomDst("GXNuclNucl4BodyMomDist", pqprC, psC, verbose) {;}
