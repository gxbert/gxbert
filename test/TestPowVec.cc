
#include "GXPow.hh"
#include "GXPowVec.hh"
#include "ApproxEqual.hh"
#include "GXInuclSpecialFunctions.hh"
#include "VecCore/Timer.h"

#include <iostream>

int main()
{
  using namespace gxbert;
  using namespace gxbert::GXInuclSpecialFunctions;

  constexpr int nvals = (2 << 3);

  GXPow            *oldpow = GXPow::GetInstance();
  GXPowVec<double> const* newpow = GXPowVec<double>::GetInstance();

  // testing low integers
  for(int i = 0; i < nvals; ++i) {
    double x = 2.*inuclRndm<double>();
    int n = round(16 * inuclRndm<double>() - 8);
    double oldres = oldpow->powN(x, n);
    double newres = newpow->PowN(x, n);
    if (i<20) {
      std::cerr<<"=== TestPowVec(i<20): "<< i <<' '<< x <<' '<< n <<' '<< oldres <<' '<< newres <<' '<<' '<< newres-oldres <<'\n';
    }
    assert( ApproxEqual(newres, oldres) );
  }

  /// vectorized version
  GXPowVec<Real_v> const* vecpow = GXPowVec<Real_v>::GetInstance();

  // quick-and-dirty benchmark
  using Real_v = VectorBackend::Double_v;
  size_t vsize = VectorSize<Real_v>();
  // nvals must be a multiple of vsize!
  assert(nvals % vsize == 0);

  unsigned int memSizeAlloc = nvals * sizeof(double);
  double *x = (double*)_mm_malloc(memSizeAlloc, 64);  // align by 64 bits
  assert(x);

  int *n = (int*)_mm_malloc(memSizeAlloc, 64);
  assert(n);

  // fill with random values (vectorized)
  for(int i = 0; i < nvals; ++i) {
    x[i] = 2. * inuclRndm<double>();
    n[i] = round(16 * inuclRndm<double>() - 8);
  }

  // benchmarkm - original
  int nReps        = 1;
  double xval, xsum = 0.;
  Timer<nanoseconds> timer;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    xsum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      xval = oldpow->powN(x[i], n[i]);
      std::cerr<<" orig bench: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< xval <<"\n";
    }
  }
  double origElapsed = timer.Elapsed();
  std::cerr<<"Original: "<< xsum <<' '<< origElapsed <<"\n";

  // benchmark - scalar
  double scx, scsum;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    scsum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      scx = newpow->PowN(x[i], n[i]);
      scsum += scx;
      std::cerr<<" scalar bench: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< scx <<'\t'<< scsum <<"\n";
    }
  }
  double scalarElapsed = timer.Elapsed();
  std::cerr<<"Scalar: "<< scsum <<' '<< scalarElapsed <<"\n";

  // benchmark - scalar v2
  double sc2x, sc2sum;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    sc2sum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      sc2x = newpow->PowN(x[i], n[i]);
      sc2sum += sc2x;
      std::cerr<<" scalar bench: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< sc2x <<'\t'<< sc2sum <<"\n";
    }
  }
  double scalar2Elapsed = timer.Elapsed();
  std::cerr<<"Scalar v2: "<< sc2sum <<' '<< scalar2Elapsed <<"\n";

  // benchmark - vectorized
  Real_v vecx, vecsum;
  Real_v const* vx = (Real_v*)x;
  Int_v const* vn = (Int_v*)n;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    scsum = 0.;
    for(size_t i = 0; i < nvals/vsize; i += vsize) {
      vecx = vecpow->PowN(vx[i], vn[i]);
      vecsum += vecx;
      std::cerr<<" vector bench: "<< i <<' '<< vx[i] <<' '<< vn[i] <<'\t'<< vecx <<'\t'<< vecsum <<"\n";
    }
  }
  double vectorElapsed = timer.Elapsed();
  std::cerr<<"Vector: "<< scsum <<' '<< vectorElapsed <<"\n";

  // cleanup
  // if(oldpow) delete oldpow;
  // if(newpow) delete newpow;
  //if(vecpow) delete vecpow;

  //=== display result
  std::cerr<<">>> PowVec tests passed.\n";
  return 0;
}
