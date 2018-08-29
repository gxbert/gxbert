//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon >= 4-body final states

#include "GXHadNucl4BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const double pqprC[2][4][4] = {
    { { 1.9439, -0.3464,  0.0271, -0.0007},
      {-4.6268,  1.1093, -0.1164,  0.0051},
      { 9.7879, -1.9313,  0.2697,  -0.015},
      {-9.6074,  1.7064, -0.3185,  0.0196} },
    { { 1.8693, -0.4996,  0.0462, -0.0013},
      {-5.5678,  1.7874, -0.1854,  0.0058},
      { 14.795,  -4.133,  0.4531, -0.0146},
      {-16.903,  3.8393, -0.4627,  0.0156} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const double psC[2][3] = {
    { 0.1491, 0.385, -0.0128 }, { 0.1802, 0.3302, -0.0094 }
  };
}

// Constructor passes arrays to templated base class

gxbert::GXHadNucl4BodyMomDst::GXHadNucl4BodyMomDst(int verbose)
  : GXVMultiBodyMomDst("GXHadNucl4BodyMomDist", pqprC, psC, verbose)
{ }
