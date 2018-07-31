
#include "GXPow.hh"
#include "GXPowVec.hh"
#include "ApproxEqual.hh"
#include "GXInuclSpecialFunctions.hh"
#include "VecCore/Timer.h"

#include <iostream>

using namespace gxbert;

constexpr size_t nReps = 1;
Timer<nanoseconds> timer;

// benchmark - original
void OriginalLoop(double* x, int* n, size_t imax) {
  GXPow *oldpow = GXPow::GetInstance();
  double xval, xsum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < imax; ++i) {
      xval = oldpow->powN(x[i], n[i]);
      xsum += xval;
      std::cerr<<" orig bench: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< xval <<' '<< xsum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<<"Original: "<< xsum <<' '<< elapsed <<"\n";
}

template <typename Real_T, typename Int_T>
void SIMDLoop(const char* testname, double const* x, int const* n, size_t nobjs)
{
  GXPowVec<Real_T> const* newpow = GXPowVec<Real_T>::GetInstance();
  std::cerr<<"VectorSizes: "<< VectorSize<Real_T>() <<" and "<< VectorSize<Int_T>() <<"\n";

  size_t imax = nobjs / VectorSize<Real_T>();
  Real_T const* xx = (Real_T const*)x;
  Int_T  const* nn = (Int_T  const*)n;
  Real_T val, sum;

  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sum = 0.;
    for(size_t i = 0; i < imax; ++i) {
      val = newpow->PowN(xx[i], nn[i]);
      sum += val;
      std::cerr<<"SIMDLoop<"<< testname <<">: "<< i <<' '<< xx[i] <<' '<< nn[i] <<'\t'<< val <<'\t'<< sum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<< testname <<": "<< sum <<' '<< elapsed <<"\n";
}


int main()
{
  using namespace gxbert;
  using namespace gxbert::GXInuclSpecialFunctions;

  constexpr int nvals = (2 << 3);

  GXPowVec<double> const* newpow = GXPowVec<double>::GetInstance();
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

  // benchmark - original
  OriginalLoop(x, n, nvals);

  // benchmark - scalar
  SIMDLoop<double,int>("Scalar", x, n, nvals);
  SIMDLoop<Real_v,gxbert::Int_v>("Vector", x, n, nvals);


  // benchmark - scalar v2
  double sc2x, sc2sum;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sc2sum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      sc2x = newpow->PowN(x[i], n[i]);
      sc2sum += sc2x;
      std::cerr<<" scalar v2: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< sc2x <<'\t'<< sc2sum <<"\n";
    }
  }
  double scalar2Elapsed = timer.Elapsed();
  std::cerr<<"Scalar v2: "<< sc2sum <<' '<< scalar2Elapsed <<"\n";

  // benchmark - vectorized
  Real_v vecx, vecsum;
  Real_v const* vx = (Real_v*)x;
  Int_v const* vn = (Int_v*)n;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    vecsum = 0.;
    for(size_t i = 0; i < nvals/vsize; i += vsize) {
      vecx = vecpow->PowN(vx[i], vn[i]);
      vecsum += vecx;
      std::cerr<<" vector bench: "<< i <<' '<< vx[i] <<' '<< vn[i] <<'\t'<< vecx <<'\t'<< vecsum <<"\n";
    }
  }
  double vectorElapsed = timer.Elapsed();
  std::cerr<<"Vector: "<< vecsum <<' '<< vectorElapsed <<"\n";

  // cleanup
  // if(oldpow) delete oldpow;
  // if(newpow) delete newpow;
  //if(vecpow) delete vecpow;

  //=== display result
  std::cerr<<">>> PowVec tests passed.\n";
  return 0;
}
