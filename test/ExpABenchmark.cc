//
// File:    ExpABenchmark.cpp
// Purpose: Benchmark for the GXPowVec::ExpA() function
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
void OriginalLoop(double* x, size_t imax) {
  GXPow *oldpow = GXPow::GetInstance();
  double xval, xsum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < imax; ++i) {
      xval = oldpow->expA(x[i]);
      xsum += xval;
      //std::cerr<<" original debug: "<< i <<' '<< x[i] <<'\t'<< xval <<' '<< xsum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<<"  Original ExpA(): "<< xsum <<' '<< elapsed/1000. <<" msec\n";
}

template <typename Real_T>
void SIMDLoop(const char* testname, const double *__restrict__ x, size_t nobjs)
{
  constexpr size_t vsize = vecCore::VectorSize<Real_T>();
  using Int_T = typename backend::VcSimdArray<vsize>::Int_v;
  GXPowVec<Real_T,Int_T> const* newpow = GXPowVec<Real_T,Int_T>::GetInstance();

  size_t imax = nobjs / VectorSize<Real_T>();
  Real_T const* xx = (Real_T const*)x;

  Real_T val, sum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < imax; ++i) {
      val = newpow->ExpA(xx[i]);
      sum += val;
      //std::cerr<<"SIMDLoop<"<< testname <<">: "<< i <<' '<< xx[i] <<'\t'<< val <<'\t'<< sum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<< testname <<" ExpA(): "<< vecCore::ReduceAdd(sum) <<' '<< elapsed/1000 <<" msec\n";
}


template <typename Real_T>
void SIMDLoopTest(const char* testname, const double *__restrict__ x, size_t nobjs)
{
  constexpr size_t vsize = vecCore::VectorSize<Real_T>();
  using Int_T = typename backend::VcSimdArray<vsize>::Int_v;
  GXPowVec<Real_T,Int_T> const* newpow = GXPowVec<Real_T,Int_T>::GetInstance();

  size_t imax = nobjs / VectorSize<Real_T>();
  Real_T const* xx = (Real_T const*)x;

  Real_T val, sum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < imax; ++i) {
      val = newpow->ExpAVec(xx[i]);  // testing ExpAVec, otherwise identical to SIMDLoop()
      sum += val;
      //std::cerr<<"SIMDLoopTest<"<< testname <<">: "<< i <<' '<< xx[i] <<'\t'<< val <<'\t'<< sum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<< testname <<" ExpAVec(): "<< vecCore::ReduceAdd(sum) <<' '<< elapsed/1000 <<" msec\n";
}


int main()
{
  constexpr int nvals = (1 << 21);

  using namespace gxbert;
  using namespace gxbert::GXInuclSpecialFunctions;

  using Real_v = VectorBackend::Double_v;
  constexpr size_t vsize = VectorSize<Real_v>();
  using Int_v = vecCore::backend::VcSimdArray<vsize>::Int_v;

  // quick-and-dirty benchmark
  // nvals must be a multiple of vsize!
  assert(nvals % vsize == 0);

  std::cerr<<"\n\n ExpA() Benchmark: nReps="<< nReps <<" and nvals="<< nvals <<"\n\n";
  unsigned int memSizeAlloc = nvals * sizeof(double);
  double *x = (double*)_mm_malloc(memSizeAlloc, 64);  // align by 64 bits
  assert(x);

  // fill with random values (vectorized)
  for(int i = 0; i < nvals; ++i) {
    //x[i] = 168. * inuclRndm<double>() + 1.;
    x[i] = 510. * inuclRndm<double>() + 1.;
  }

  // benchmark - original
  OriginalLoop(x, nvals);

  // benchmark - scalar
  SIMDLoop<double>("    Scalar", x, nvals);

  // benchmark - vector
  SIMDLoop<Real_v>("    Vector", x, nvals);

  // benchmark - scalar type + std:: function

  // benchmark - scalar
  SIMDLoopTest<double>(" Scalar", x, nvals);

  // benchmark - vector
  SIMDLoopTest<Real_v>(" Vector", x, nvals);

  // benchmark - scalar type + std:: function
  GXPowVec<double, int>    const* newpow = GXPowVec<double,int>::GetInstance();
  double sc2x, sc2sum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < nvals; ++i) {
      sc2x = newpow->StdExp(x[i]);
      sc2sum += sc2x;
      //std::cerr<<" scalarStd debug: "<< i <<' '<< a[i] <<' '<< x[i] <<'\t'<< sc2x <<'\t'<< sc2sum <<"\n";
    }
  }
  double scalar2Elapsed = timer.Elapsed();
  std::cerr<<" ScalarStd ExpA(): "<< sc2sum <<' '<< scalar2Elapsed/1000. <<" msec\n";

  // benchmark - vectorized type + std:: function
  GXPowVec<Real_v, Int_v> const* vecpow = GXPowVec<Real_v,Int_v>::GetInstance();
  Real_v vecx, vecsum = 0.;
  Real_v const* vx = (Real_v*)x;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < nvals/vsize; ++i) {
      vecx = vecpow->StdExp(vx[i]);
      vecsum += vecx;
      //std::cerr<<" vectorStd debug: "<< i <<' '<< va[i] <<' '<< vx[i] <<'\t'<< vecx <<'\t'<< vecsum <<"\n";
    }
  }
  double vectorElapsed = timer.Elapsed();
  std::cerr<<" VectorStd ExpA(): "<< vecCore::ReduceAdd(vecsum) <<' '<< vectorElapsed/1000. <<" msec\n";

  return 0;
}
