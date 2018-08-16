//
// File:    TestPowVec.cpp
// Purpose: Unit tests for the GXPowVec functions
//
// 2018-07-20 Guilherme Lima  - Created
//

//.. ensure asserts are compiled in
#undef NDEBUG

#include "GXPow.hh"
#include "GXPowVec.hh"
#include "ApproxEqual.hh"
#include "GXInuclSpecialFunctions.hh"
#include "VecCore/Timer.h"

#include <iostream>

#define DEBUG 1

int main()
{
  using namespace gxbert;
  using namespace gxbert::GXInuclSpecialFunctions;
  using Int_v = vecCore::backend::VcSimdArray<2>::Int_v;

  constexpr int nvals = (1 << 11);
  constexpr int maxZ = 512;

  GXPow const* oldpow = GXPow::GetInstance();
  GXPowVec<double, int> const* newpow = GXPowVec<double,int>::GetInstance();

  /// check vectorized version can be instantiated
  GXPowVec<Real_v, Int_v> const* vecpow = GXPowVec<Real_v, Int_v>::GetInstance();
  (void)vecpow;  // avoid compilation warning due to unused variable

  //===========
  std::cerr <<"--- testing PowN(x,N) with low integers...\n";
  for(int i = 0; i < nvals; ++i) {
    double x = 40.*inuclRndm<double>() - 20.;
    int n = round(16 * inuclRndm<double>() - 8);
    double stdval = std::pow(x, n);
    double oldval = oldpow->powN(x, n);
    double newval = newpow->PowN(x, n);

#ifdef DEBUG
    //if (i<20) {
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"*** Problem: "<< i <<' '<< x <<' '<< n <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }

  //==== testing PowN(x,N) for larger ranges of integers
  std::cerr <<"--- testing PowN(x,N) for larger ranges of integers...\n";
  for(int i = 0; i < nvals; ++i) {
    double x = 2.*inuclRndm<double>();
    int n = round(30 * inuclRndm<double>() - 15);
    double stdval = std::pow(x, n);
    double oldval = oldpow->powN(x, n);
    double newval = newpow->PowN(x, n);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"*** Problem: "<< i <<' '<< x <<' '<< n <<' '<< stdval <<" "<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }

  //==== testing Z13()
  std::cerr <<"--- testing Z13(n)...\n";
  constexpr double oneThird = 1./3.;
  for(int n = 1; n < maxZ; ++n) {
    double stdval = std::pow(n, oneThird);
    double oldval = oldpow->Z13(n);
    double newval = newpow->Z13(n);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: Z13() problem: n="<< n <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    assert( ApproxEqual(newval, stdval) );
  }

  //==== testing A13()
  std::cerr <<"--- testing A13(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double a = 512. * double(i) / double(nvals);
    double stdval = std::pow(a, oneThird);
    double oldval = oldpow->A13(a);
    double newval = newpow->A13(a);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: A13() problem: a="<< a <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }


  //==== testing Z23()
  std::cerr <<"--- testing Z23(n)...\n";
  constexpr double twoThirds = 2./3.;
  for(int n = 1; n < maxZ; ++n) {
    double stdval = std::pow(n, twoThirds);
    double oldval = oldpow->Z23(n);
    double newval = newpow->Z23(n);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: Z23() problem: n="<< n <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    assert( ApproxEqual(newval, stdval) );
  }

  //==== testing A23()
  std::cerr <<"--- testing A23(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double a = 512. * double(i) / double(nvals);
    double stdval = std::pow(a, twoThirds);
    double oldval = oldpow->A23(a);
    double newval = newpow->A23(a);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: A23() problem: a="<< a <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }


  //==== testing LogZ()
  std::cerr <<"--- testing LogZ(N)...\n";
  for(int n = 1; n < maxZ; ++n) {
    double stdval = std::log(n);
    double oldval = oldpow->logZ(n);
    double newval = newpow->LogZ(n);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: LogZ() problem: n="<< n <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    assert( ApproxEqual(newval, stdval) );
  }

  //==== testing LogX()
  std::cerr <<"--- testing LogX(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double x = 512. * double(i) / double(nvals);
    double stdval = std::log(x);
    double oldval = oldpow->logX(x);
    double newval = newpow->LogX(x);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: LogX() problem: x="<< x <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }

  //==== testing LogA()
  std::cerr <<"--- testing LogA(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double x = 512. * double(i) / double(nvals);
    double stdval = std::log(x);
    double oldval = oldpow->logA(x);
    double newval = newpow->LogA(x);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: LogA() problem: x="<< x <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }

  //==== testing Log10Z()
  std::cerr <<"--- testing Log10Z(N)...\n";
  const double toLogBase10 = 1.0/std::log(10.);
  for(int n = 1; n < maxZ; ++n) {
    double stdval = std::log(n) * toLogBase10;
    double oldval = oldpow->log10Z(n);
    double newval = newpow->Log10Z(n);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: Log10Z() problem: n="<< n <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error"<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    assert( ApproxEqual(newval, stdval) );
  }


  //==== testing Log10A()
  std::cerr <<"--- testing Log10A(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double x = 512. * double(i) / double(nvals);
    double stdval = std::log(x) * toLogBase10;
    double oldval = oldpow->log10A(x);
    double newval = newpow->Log10A(x);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: Log10A() problem: x="<< x <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }


  //==== testing ExpA()
  std::cerr <<"--- testing ExpA(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double x = 512. * double(i) / double(nvals);
    double stdval = std::exp(x);
    double oldval = oldpow->expA(x);
    double newval = newpow->ExpA(x);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: ExpA() problem: i="<< i <<" a="<< x <<' '<< stdval <<' '<< oldval <<' '<< newval <<' '<< newval-oldval <<" frac.error="<< (newval-stdval)/stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }


  //==== testing PowZ()
  std::cerr <<"--- testing PowZ(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    int z = round(510. * inuclRndm<double>() + 2.);
    double x = 10.*inuclRndm<double>() - 20.;
    double stdval = std::pow(z,x);
    double oldval = oldpow->powZ(z,x);
    double newval = newpow->PowZ(z,x);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: PowZ() problem: z="<< z <<", x="<< x <<", vals: "<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }


  //==== testing PowA()
  std::cerr <<"--- testing PowA(x)...\n";
  for(int i = 1; i < nvals; ++i) {
    double a = 512. * double(i) / double(nvals);
    double x = 10.*inuclRndm<double>() - 20.;
    double stdval = std::pow(a,x);
    double oldval = oldpow->powA(a,x);
    double newval = newpow->PowA(a,x);

#ifdef DEBUG
    if (std::fabs(stdval-newval) > 1.0e-6) {
      std::cerr<<"=== TestPowVec: PowA() problem: a="<< a <<", x="<< x <<", vals: "<< stdval <<' '<< oldval <<' '<< newval <<' '<<' '<< newval-oldval <<" error="<< newval-stdval <<'\n';
    }
#endif

    assert( ApproxEqual(newval, oldval) );
    //assert( ApproxEqual(newval, stdval) );
  }


  //=== display result
  std::cerr<<">>> PowVec tests passed.\n";
  return 0;
}
