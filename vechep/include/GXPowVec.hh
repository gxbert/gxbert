//
// @File: GXPowVec.hh
// -------------------------------------------------------------------
//
// Class GXPowVec
//
// Class description:
//
// Utility singleton class for the fast computation of log and pow
// functions. Integer argument should in the interval 0-512, no
// check is performed inside these methods for performance reasons.
// For factorial integer argument should be in the interval 0-170
// Computations with double arguments are fast for the interval
// 0.002-511.5 for all functions except exponent, which is computed 
// for the interval 0-84.4, standard library is used in the opposite case
//
// 2009-05-23 Vladimir Ivanchenko  Created
// 2018-07-17 Guilherme Lima       SIMD-vectorization, with name change to avoid name clash
// -------------------------------------------------------------------

#ifndef GXPowVec_hh
#define GXPowVec_hh 1

//#include "VecCore/VecCore"
#include "VecHepDefs.hh"
#include "GXPow.hh"
#include "GXLog.hh"
#include "GXExp.hh"
#include <vector>

namespace gxbert {

inline namespace
GX_IMPL_NAMESPACE {

#define GXBERT_UNARY_FUNCTION(F, f, T1)       \
VECCORE_FORCE_INLINE VECCORE_ATT_HOST_DEVICE  \
T F(const T1 &x) const                        \
{                                             \
  T ret;                                      \
  auto ppow = GXPow::GetInstance();           \
  for(size_t i = 0; i < fTsize; ++i)          \
    Set(ret, i, ppow->f(Get(x,i)));           \
  return ret;                                 \
}

#define GXBERT_BINARY_FUNCTION(F, f, T1, T2)  \
VECCORE_FORCE_INLINE VECCORE_ATT_HOST_DEVICE  \
T F(const T1 &x, const T2 &y) const	      \
{                                             \
  T ret;                                      \
  auto ppow = GXPow::GetInstance();           \
  for(size_t i = 0; i < fTsize; ++i)          \
    Set(ret, i, ppow->f(Get(x,i), Get(y,i))); \
  return ret;                                 \
}

template <typename T, typename Int_T>
class GXPowVec
{
  using Bool_v = typename vecCore::Mask_v<T>;
  static constexpr size_t fTsize = vecCore::VectorSize<T>();

public:

  static GXPowVec* GetInstance()
  {
    return fpInstance;
  }

  // Fast computation of Z^1/3
  //
  //VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T Z13(Int_T Z) const;
  GXBERT_UNARY_FUNCTION(Z13, Z13, Int_T)
  GXBERT_UNARY_FUNCTION(A13, A13, T)

  //VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T A13(T A) const;

  /// Fast computation of Z^2/3
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T Z23(Int_T Z) const
  {
    T x = Z13(Z);
    return x*x;
  }

  /// Fast computation of A^2/3
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T A23(T A) const
  {
    T x = A13(A);
    return x*x;
  }

  // Fast computation of log(Z)
  //
  // VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logZ(Int_T Z) const;
  // VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logX(T x) const;
  GXBERT_UNARY_FUNCTION(LogZ, logZ, Int_T)
  GXBERT_UNARY_FUNCTION(LogX, logX, T)

  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T LogA(T A) const
  {
    // return 1.0 <= A ? LogBase(A) : -LogBase(1./A);

    // modified version: calls LogBase() only once
    Bool_v aboveOne = A >= 1.0;
    T easierA = Blend( aboveOne, A, T(1.0)/A);
    T result = LogBase( easierA );
    MaskedAssign( result, !aboveOne, -result );
    return result;
  }

  // Fast computation of log10(Z)
  //
  // VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T log10Z(Int_T Z) const;
  // VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T log10A(T A) const;
  GXBERT_UNARY_FUNCTION(Log10Z, log10Z, Int_T)
  GXBERT_UNARY_FUNCTION(Log10A, log10A, T)

  // Fast computation of exp(X)
  //
  //VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T expA(T A) const;
  GXBERT_UNARY_FUNCTION(ExpA, expA, T)

  // Fast computation of pow(Z,X)
  //
  //VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T powZ(Int_T Z, T y) const;
  GXBERT_BINARY_FUNCTION(PowZ, powZ, Int_T, T)

  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T PowA(T A, T y) const
  {
    const T zero(0.0);
    return vecCore::Blend(A == zero ? zero : expA(y*logX(A)));
  }

  //template <typename T1>
  //VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T PowN(T x, T1 n) const;
  GXBERT_BINARY_FUNCTION(PowN, powN, T, Int_T)

  // Fast factorial
  //
  // VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T factorial(Int_T Z) const;
  // VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logfactorial(Int_T Z) const;
  GXBERT_UNARY_FUNCTION(Factorial,    factorial,    Int_T)
  GXBERT_UNARY_FUNCTION(LogFactorial, logfactorial, Int_T)

private:

  GXPowVec() {}

  ~GXPowVec()
  {
    if (GXPowVec<T,Int_T>::fpInstance) {
      delete GXPowVec<T, Int_T>::fpInstance;
      GXPowVec<T,Int_T>::fpInstance = NULL;
    }
  }

  //VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T LogBase(T x) const;
  GXBERT_UNARY_FUNCTION(LogBase, logBase, T)

  template <size_t IMAX>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T PowLowNLoop(const T xx, Index_v<T> nn) const;

  static GXPowVec<T,Int_T>* fpInstance;
};

template <typename T, typename Int_T>
GXPowVec<T, Int_T>* GXPowVec<T,Int_T>::fpInstance = new GXPowVec<T,Int_T>();

// -------------------------------------------------------------------

// returns pow(x,N) only for lanes with positive nn <= IMAX, and 1.0 otherwise
template <typename T, typename Int_T>
template <size_t IMAX> // default should be IMAX=8, but default parameter cannot be used here
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T,Int_T>::
PowLowNLoop(const T xx, gxbert::Index_v<T> nn) const
{
  size_t imax = ReduceMax(nn);
  if (imax > IMAX) imax = IMAX;
  Bool_v lowN = (nn <= Index_v<T>(IMAX));

  T res(1.0);
  Index_v<T> one(1);
  // std::cerr<<"PowLowNLoop: xx="<< xx <<" nn="<< nn <<" res="<< res <<"\n";
  for(size_t i = 0; i < imax; ++i) {
    vecCore::MaskedAssign( res, lowN && nn > 0, res * xx); // only as needed in each lane
    vecCore::MaskedAssign( nn, lowN && nn > 0, nn - one);    // keep track # times / lane
    // std::cerr<<"PowLowNLoop: i="<< i <<" nn="<< nn <<", res="<< res <<"\n";
  }
  // std::cerr<<"PowLowNLoop: return res="<< res <<"\n";
  return res;
}

/*// IT represents an integer type of same lenght in bits as type T
template <typename T, typename Int_T>
template <typename IT>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T,Int_T>::
PowN(T x, IT n) const
{
  using vecCore::Convert;

  // std::cerr<<"-- PowN() called with x="<< x <<" and n="<< n <<")\n";
  // std::cerr<<"-- PowN(): MaskedAssign("<< xx <<", "<< n <<"<0, T(1.)/"<< x <<");\n";
  vecCore::MaskedAssign( x, vecCore::Convert<T,IT>(n) < T(0.), T(1.)/x );
  vecCore::MaskedAssign( n, n < 0, -n );
  Index_v<T> nn = Abs(n);

  // use PowLowNLoop() whenever possible
  // std::cerr<<"-- PowN(): calling powLowNLoop<8>("<< xx <<", "<< nn <<"\n";
  T result = PowLowNLoop<8>(x, nn);
  // std::cerr<<"PowVec::PowN(): xx="<< x <<" "<< xx <<", nn="<< n <<" "<< nn <<", result="<< result <<"\n";

  Bool_v bigN = (Convert<T,IT>(n) > 8);
  if (!vecCore::MaskEmpty(bigN)) {
    for(size_t i = 0; i < fTsize; ++i) {
      if( Get(bigN,i) ) Set(result, i, std::pow(Get(x,i), Get(n,i)));
    }
  }
  return result;
  }*/

// -------------------------------------------------------------------
#undef GXBERT_UNARY_FUNCTION
#undef GXBERT_BINARY_FUNCTION

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace
#endif
