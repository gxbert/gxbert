//
// File:    Log10ZBenchmark.cpp
// Purpose: Benchmark for the GXPowVec::Log10Z() function
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
void OriginalLoop(int* z, size_t imax) {
  GXPow *oldpow = GXPow::GetInstance();
  double xval, xsum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    xsum = 0.;
    for(size_t i = 0; i < imax; ++i) {
      xval = oldpow->log10Z(z[i]);
      xsum += xval;
      //std::cerr<<" original debug: "<< i <<' '<< z[i] <<'\t'<< xval <<' '<< xsum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<<"  Original Log10Z(): "<< xsum <<' '<< elapsed/1000. <<" msec\n";
}

template <typename Real_T>
void SIMDLoop(const char* testname, const int *__restrict__ z, size_t nobjs)
{
  constexpr size_t vsize = vecCore::VectorSize<Real_T>();
  using Int_T = typename backend::VcSimdArray<vsize>::Int_v;
  GXPowVec<Real_T,Int_T> const* newpow = GXPowVec<Real_T,Int_T>::GetInstance();

  size_t imax = nobjs / VectorSize<Real_T>();
  Int_T const* zz = (Int_T const*)z;

  Real_T val, sum;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sum = 0.;
    for(size_t i = 0; i < imax; ++i) {
      val = newpow->Log10Z(zz[i]);
      sum += val;
      //std::cerr<<"SIMDLoop<"<< testname <<">: "<< i <<' '<< zz[i] <<'\t'<< val <<'\t'<< sum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<< testname <<" Log10Z(): "<< vecCore::ReduceAdd(sum) <<' '<< elapsed/1000 <<" msec\n";
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

  std::cerr<<"\n\n Log10Z() Benchmark: nReps="<< nReps <<" and nvals="<< nvals <<"\n\n";
  unsigned int memSizeAlloc = nvals * sizeof(int);
  int *z = (int*)_mm_malloc(memSizeAlloc, 64);  // align by 64 bits
  assert(z);

  // fill with random values (vectorized)
  for(int i = 0; i < nvals; ++i) {
    z[i] = round(169. * inuclRndm<double>() + 1.);
  }

  // benchmark - original
  OriginalLoop(z, nvals);

  // benchmark - scalar
  SIMDLoop<double>("    Scalar", z, nvals);

  // benchmark - vector
  SIMDLoop<Real_v>("    Vector", z, nvals);

  // benchmark - scalar type + std:: function
  GXPowVec<double, int>    const* newpow = GXPowVec<double,int>::GetInstance();
  double sc2x, sc2sum;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sc2sum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      sc2x = newpow->StdLog10(z[i]);
      sc2sum += sc2x;
      //std::cerr<<" scalarStd debug: "<< i <<' '<< z[i] <<'\t'<< sc2x <<'\t'<< sc2sum <<"\n";
    }
  }
  double scalar2Elapsed = timer.Elapsed();
  std::cerr<<" ScalarStd Log10Z(): "<< sc2sum <<' '<< scalar2Elapsed/1000. <<" msec\n";

  // benchmark - vectorized type + std:: function
  GXPowVec<Real_v, Int_v> const* vecpow = GXPowVec<Real_v,Int_v>::GetInstance();
  Real_v vecx, vecsum;
  Int_v const* vz = (Int_v*)z;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    vecsum = 0.;
    for(size_t i = 0; i < nvals/vsize; ++i) {
      vecx = vecpow->StdLog10(vz[i]);
      vecsum += vecx;
      //std::cerr<<" vectorStd debug: "<< i <<' '<< vz[i] <<'\t'<< vecx <<'\t'<< vecsum <<"\n";
    }
  }
  double vectorElapsed = timer.Elapsed();
  std::cerr<<" VectorStd Log10Z(): "<< vecCore::ReduceAdd(vecsum) <<' '<< vectorElapsed/1000. <<" msec\n";

  return 0;
}
