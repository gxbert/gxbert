
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
void OriginalLoop(double* x, int* n, size_t imax) {
  GXPow *oldpow = GXPow::GetInstance();
  double xval, xsum = 0.;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    for(size_t i = 0; i < imax; ++i) {
      xval = oldpow->powN(x[i], n[i]);
      xsum += xval;
      //std::cerr<<" orig bench: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< xval <<' '<< xsum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<<"  Original: "<< xsum <<' '<< elapsed/1000. <<" msec\n";
}

template <typename Real_T, typename Int_T>
void SIMDLoop(const char* testname, double const* x, int const* n, size_t nobjs)
{
  GXPowVec<Real_T,Int_T> const* newpow = GXPowVec<Real_T,Int_T>::GetInstance();

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
      //std::cerr<<"SIMDLoop<"<< testname <<">: "<< i <<' '<< xx[i] <<' '<< nn[i] <<'\t'<< val <<'\t'<< sum <<"\n";
    }
  }
  double elapsed = timer.Elapsed();
  std::cerr<< testname <<": "<< vecCore::ReduceAdd(sum) <<' '<< elapsed/1000 <<" msec\n";
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

  std::cerr<<"\n\n PowVec Benchmark: nReps="<< nReps <<" and nvals="<< nvals <<"\n\n";
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

  // benchmark - scalar v0
  GXPowVec<double, int>    const* newpow = GXPowVec<double,int>::GetInstance();
  double sc2x, sc2sum;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    sc2sum = 0.;
    for(size_t i = 0; i < nvals; ++i) {
      sc2x = newpow->PowN(x[i], n[i]);
      sc2sum += sc2x;
      //std::cerr<<" scalar v0 bench: "<< i <<' '<< x[i] <<' '<< n[i] <<'\t'<< sc2x <<'\t'<< sc2sum <<"\n";
    }
  }
  double scalar2Elapsed = timer.Elapsed();
  std::cerr<<" Scalar v0: "<< sc2sum <<' '<< scalar2Elapsed/1000. <<" msec\n";

  // benchmark - vectorized
  GXPowVec<Real_v, Int_v> const* vecpow = GXPowVec<Real_v,Int_v>::GetInstance();
  Real_v vecx, vecsum;
  Real_v const* vx = (Real_v*)x;
  Int_v const* vn = (Int_v*)n;
  timer.Start();
  for(size_t irep = 0; irep < nReps; ++irep) {
    vecsum = 0.;
    for(size_t i = 0; i < nvals/vsize; ++i) {
      vecx = vecpow->PowN(vx[i], vn[i]);
      vecsum += vecx;
      //std::cerr<<" vector v0 bench: "<< i <<' '<< vx[i] <<' '<< vn[i] <<'\t'<< vecx <<'\t'<< vecsum <<"\n";
    }
  }
  double vectorElapsed = timer.Elapsed();
  std::cerr<<" Vector v0: "<< vecCore::ReduceAdd(vecsum) <<' '<< vectorElapsed/1000. <<" msec\n";

  // cleanup
  // if(oldpow) delete oldpow;
  // if(newpow) delete newpow;
  //if(vecpow) delete vecpow;

  // benchmark - scalar
  SIMDLoop<double, int>("    Scalar", x, n, nvals);

  // benchmark - vector
  SIMDLoop<Real_v, Int_v>("    Vector", x, n, nvals);

  return 0;
}