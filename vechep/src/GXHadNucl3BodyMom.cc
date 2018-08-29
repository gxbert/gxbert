//
// Description: class containing parametrized momentum distributions
//              in the CM for hadron/nucleon 3-body final states
//
#include "GXHadNucl3BodyMomDst.hh"

namespace {
  // Powers of Ekin^0..3, blocks of S^0..3 for PQ,PR: outgoing N; outgoing h,K,Y
  static const double pqprC[2][4][4] = {
    { { 0.6305,  2.1801, -1.2886,  0.2091},
      {-3.7333,  1.5163,  -2.457,  0.5228},
      { 13.464,  -16.38,  15.129, -2.8687},
      {-18.594,  27.944, -23.295,  4.2688} },
    { { 0.9336,  1.7811, -1.5264,  0.2713},
      {-1.8181, -8.2927,  6.8433, -1.1944},
      { 5.5157,  20.607, -16.067,  2.7495},
      {-8.5216, -20.827,  16.845, -2.9045} }
  };

  // Powers of Ekin^0..2 for PS: outgoing N; outgoing h,K,Y
  static const double psC[2][3] = {
    { 0.0929,  0.5389, -0.0545 }, { 0.1303,  0.4071, -0.0288 }
  };
}

// Constructor passes arrays to templated base class

gxbert::GXHadNucl3BodyMomDst::GXHadNucl3BodyMomDst(int verbose)
  : GXVMultiBodyMomDst("GXHadNucl3BodyMomDist", pqprC, psC, verbose)
{ }
