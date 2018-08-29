//
// Description: class containing parametrized momentum distributions
//              in the CM for nucleon/nucleon 3-body final states
//
#include "GXNuclNucl3BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const double pqprC[2][4][4] = {
    { { 0.5028,  0.9348, -0.0967,  -0.025},
      { 3.1442,  -10.59,  4.7335, -0.6248},
      {-7.8172,  29.227, -14.298,  2.0282},
      { 8.1667,  -34.55,  17.685, -2.5895} },
    { { 1.1965,   0.287, -0.2449,  0.0373},
      {-0.8289, -4.9065,  2.9191,  -0.422},
      { 1.0426,  16.264, -9.5776,  1.3883},
      { -1.909, -19.904,  11.938, -1.7476} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const double psC[2][3] = {
    { 0.1451,  0.4652,  -0.033 }, { 0.1538,  0.2744, -0.0146 }
  };
}

// Constructor passes arrays to templated base class

gxbert::GXNuclNucl3BodyMomDst::GXNuclNucl3BodyMomDst(int verbose)
  : GXVMultiBodyMomDst("GXNuclNucl3BodyMomDist", pqprC, psC, verbose) {;}
