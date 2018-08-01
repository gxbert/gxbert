//
// File:    PowABenchmark.cpp
// Purpose: Benchmark for the GXPowVec::PowA() function
//
// 2018-08-01 Guilherme Lima  - Created
//

#include "GXPow.hh"
#include "GXPowVec.hh"
#include "ApproxEqual.hh"
#include "GXInuclSpecialFunctions.hh"
#include "VecCore/Timer.h"
#include "VecCore/VecCore"

#include <iostream>

using namespace gxbert;

constexpr size_t nReps = 100;
Timer<microseconds> timer;

// benchmark - original
void OriginalLoop(double* a, double* x, size_t imax) {
  GXPow *oldpow = GXPow::GetInstance();
  double xval, xsum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < imax; ++i) {
      xval = oldpow->powA(a[i], x[i]);
      xsum += xval;
      //std::cerr<<" original debug: "<< i <<' '<< a[i] <<' '<< x[i] <<'\t'<< xval <<' '<< xsum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<<"  Original PowA(): "<< xsum <<' '<< elapsed/1000. <<" msec\n";
}

template <typename Real_T>
void SIMDLoop(const char* testname, const double *__restrict__ a, const double *__restrict__ x, size_t nobjs)
{
  constexpr size_t vsize = vecCore::VectorSize<Real_T>();
  using Int_T = typename backend::VcSimdArray<vsize>::Int_v;
  GXPowVec<Real_T,Int_T> const* newpow = GXPowVec<Real_T,Int_T>::GetInstance();

  size_t imax = nobjs / VectorSize<Real_T>();
  Real_T const* aa = (Real_T const*)a;
  Real_T const* xx = (Real_T const*)x;

  Real_T val, sum;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sum = 0.;
    for(size_t i = 0; i < imax; ++i) {
      val = newpow->PowA(aa[i], xx[i]);
      sum += val;
      //std::cerr<<"SIMDLoop<"<< testname <<">: "<< i <<' '<< xx[i] <<' '<< aa[i] <<'\t'<< val <<'\t'<< sum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<< testname <<" PowA(): "<< vecCore::ReduceAdd(sum) <<' '<< elapsed/1000 <<" msec\n";
}


int main()
{
  constexpr int nvals = (2 << 20);

  using namespace gxbert;
  using namespace gxbert::GXInuclSpecialFunctions;

  using Real_v = VectorBackend::Double_v;
  constexpr size_t vsize = VectorSize<Real_v>();
  using Int_v = vecCore::backend::VcSimdArray<vsize>::Int_v;

  // quick-and-dirty benchmark
  // nvals must be a multiple of vsize!
  assert(nvals % vsize == 0);

  std::cerr<<"\n\n PowA() Benchmark: nReps="<< nReps <<" and nvals="<< nvals <<"\n\n";
  unsigned int memSizeAlloc = nvals * sizeof(double);
  double *x = (double*)_mm_malloc(memSizeAlloc, 64);  // align by 64 bits
  assert(x);

  double *a = (double*)_mm_malloc(memSizeAlloc, 64);
  assert(a);

  // fill with random values (vectorized)
  for(int i = 0; i < nvals; ++i) {
    a[i] = round(510. * inuclRndm<double>() + 1.);
    x[i] = 10. * inuclRndm<double>() - 5.;
  }

  // benchmark - original
  OriginalLoop(a, x, nvals);

  // benchmark - scalar v0
  GXPowVec<double, int>    const* newpow = GXPowVec<double,int>::GetInstance();
  double sc2x, sc2sum;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sc2sum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      sc2x = newpow->StdPowA(a[i], x[i]);
      sc2sum += sc2x;
      //std::cerr<<" scalar v0 debug: "<< i <<' '<< a[i] <<' '<< x[i] <<'\t'<< sc2x <<'\t'<< sc2sum <<"\n";
    }
  }
  double scalar2Elapsed = timer.Elapsed();
  std::cerr<<" Scalar v0 PowA(): "<< sc2sum <<' '<< scalar2Elapsed/1000. <<" msec\n";

  // benchmark - vectorized
  GXPowVec<Real_v, Int_v> const* vecpow = GXPowVec<Real_v,Int_v>::GetInstance();
  Real_v vecx, vecsum;
  Int_v const* va = (Int_v*)a;
  Real_v const* vx = (Real_v*)x;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    vecsum = 0.;
    for(size_t i = 0; i < nvals/vsize; ++i) {
      vecx = vecpow->StdPowA(va[i], vx[i]);
      vecsum += vecx;
      //std::cerr<<" vector v0 debug: "<< i <<' '<< va[i] <<' '<< vx[i] <<'\t'<< vecx <<'\t'<< vecsum <<"\n";
    }
  }
  double vectorElapsed = timer.Elapsed();
  std::cerr<<" Vector v0 PowA(): "<< vecCore::ReduceAdd(vecsum) <<' '<< vectorElapsed/1000. <<" msec\n";

  // benchmark - scalar
  SIMDLoop<double>("    Scalar", a, x, nvals);

  // benchmark - vector
  SIMDLoop<Real_v>("    Vector", a, x, nvals);

  return 0;
}
