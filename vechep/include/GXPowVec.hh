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

//#include "globals.hh"
#include "GXLog.hh"
#include "GXExp.hh"
#include <vector>

namespace gxbert {

inline namespace
GX_IMPL_NAMESPACE {

template <typename T>
class GXPowVec
{
  using Int_v  = vecCore::Index_v<T>;
  using Bool_v = vecCore::Mask_v<T>;

public:

  ~GXPowVec();

  static GXPowVec* GetInstance()
  {
    if (fpInstance == 0) {
      fpInstance = new GXPowVec<T>();
    }
    return fpInstance;
  }

  // Fast computation of Z^1/3
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T Z13(Int_v Z) const;

  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T A13(T A) const;

  // Fast computation of Z^2/3
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T Z23(Int_v Z) const;

  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T A23(T A) const;

  // Fast computation of log(Z)
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logZ(Int_v Z) const;
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logA(T A) const;
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logX(T x) const;

  // Fast computation of log10(Z)
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T log10Z(Int_v Z) const;
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T log10A(T A) const;

  // Fast computation of exp(X)
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T expA(T A) const;

  // Fast computation of pow(Z,X)
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T powZ(Int_v Z, T y) const;
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T powA(T A, T y) const;
  T powN(T x, Int_v n) const;

  // Fast factorial
  //
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T factorial(Int_v Z) const;
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logfactorial(Int_v Z) const;

private:

  GXPowVec();

  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T logBase(T x) const;

  template <size_t IMAX>
  VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T powLowNLoop(const T xx, Int_v nn) const;

  static GXPowVec<T>* fpInstance;

  const T onethird;
  const Int_v    max2;

  T maxA;
  T maxA2;
  T maxAexp;

  std::vector<T> ener;
  std::vector<T> logen;
  std::vector<T> pz13;
  std::vector<T> lz;
  std::vector<T> lz2;
  std::vector<T> fexp;
  std::vector<T> fact;
  std::vector<T> logfact;
};

template <typename T>
GXPowVec<T>* GXPowVec<T>::fpInstance = new GXPowVec<T>();

// -------------------------------------------------------------------

template <typename T> GXPowVec<T>::
GXPowVec()
  : onethird(1.0/3.0), max2(5)
{  
  const G4int maxZ = 512; 
  const G4int maxZfact = 170; 

  maxA    = -0.6 + maxZ;
  maxA2   = 1.25 + max2*0.2;
  maxAexp = -0.76+ maxZfact*0.5;

  ener.resize(max2+1,1.0);
  logen.resize(max2+1,0.0);
  lz2.resize(max2+1,0.0);
  pz13.resize(maxZ,0.0);
  lz.resize(maxZ,0.0);
  fexp.resize(maxZfact,0.0);
  fact.resize(maxZfact,0.0);
  logfact.resize(maxZ,0.0);

  G4double f = 1.0;
  G4double logf = 0.0;
  fact[0] = 1.0;
  fexp[0] = 1.0;

  for(G4int i=1; i<=max2; ++i)
  {
    ener[i] = powN(500.,i); 
    logen[i]= GXLog(ener[i]); 
    lz2[i]  = GXLog(1.0 + i*0.2);
  }

  for(G4int i=1; i<maxZ; ++i)
  {
    G4double x  = G4double(i);
    pz13[i] = std::pow(x,onethird);
    lz[i]   = GXLog(x);
    if(i < maxZfact)
    { 
      f *= x; 
      fact[i] = f;
      fexp[i] = GXExp(0.5*i);
    }
    logf += lz[i];
    logfact[i] = logf;
  }
}

template <typename T>
GXPowVec<T>::~GXPowVec()
{
  if (GXPowVec<T>::fpInstance) {
    delete GXPowVec<T>::fpInstance;
    GXPowVec<T>::fpInstance = NULL;
  }
}
  
template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
Z13(Int_v Z) const
{
  return pz13[Z];
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
A13(T A) const
{
  T res = 0.0;
  if(A > 0.0) 
  {
    T a = (1.0 <= A) ? A : 1.0/A;
    if(1.0 > A) { a = 1.0/A; }
    if(a <= maxA)
    {
      Int_v i = Int_v(a + 0.5);
      T x = (a/T(i) - 1.0)*onethird;
      res = pz13[i]*(1.0 + x - x*x*(1.0 - 1.66666666*x));
      if(1.0 > A) { res = 1.0/res; }
    }
    else
    {
      res = std::pow(A, onethird); 
    }
  }
  return res;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
Z23(Int_v Z) const
{
  T x = Z13(Z);
  return x*x;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
A23(T A) const
{
  T x = A13(A);
  return x*x;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
logZ(Int_v Z) const
{
  return lz[Z];
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
logBase(T a) const
{
  T res;
  if(a <= maxA2) 
  {
    Int_v i = Int_v(max2*(a - 1) + 0.5);
    if(i > max2) { i = max2; }
    T x = a/(T(i)/max2 + 1) - 1;
    res = lz2[i] + x*(1.0 - (0.5 - onethird*x)*x);
  }
  else if(a <= maxA)
  {
    Int_v i = Int_v(a + 0.5);
    T x = a/T(i) - 1;
    res = lz[i] + x*(1.0 - (0.5 - onethird*x)*x);
  }
  else
  {
    res = GXLog(a);
  }
  return res;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
logA(T A) const
{
  return (1.0 <= A ? logBase(A) : -logBase(1./A));
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
logX(T x) const
{
  T res = 0.0;
  T a = (1.0 <= x) ? x : 1.0/x;

  if(a <= maxA) 
  {
    res = logBase(a);
  }
  else if(a <= ener[2])
  {
    res = logen[1] + logBase(a/ener[1]);
  }
  else if(a <= ener[3])
  {
    res = logen[2] + logBase(a/ener[2]);
  }
  else
  {
    res = GXLog(a);
  }

  if(1.0 > x) { res = -res; }
  return res;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
log10Z(Int_v Z) const
{
  return lz[Z]/lz[10];
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
log10A(T A) const
{
  return logX(A)/lz[10];
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
expA(T A) const
{
  T res;
  T a = (0.0 <= A) ? A : -A;

  if(a <= maxAexp)
  {
    Int_v i = Round(2.*a);
    T x = a - 0.5 * T(i);
    res = fexp[i]*(1.0 + x*(1.0 + 0.5*(1.0 + onethird*x)*x));
  }
  else
  {
    res = GXExp(a);
  }
  if(0.0 > A) { res = 1.0/res; }
  return res;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
powZ(Int_v Z, T y) const
{
  return expA(y*lz[Z]);
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
powA(T A, T y) const
{
  return (0.0 == A ? 0.0 : expA(y*logX(A)));
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
factorial(Int_v Z) const
{
  return fact[Z];
}

template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
logfactorial(Int_v Z) const
{
  return logfact[Z];
}

// -------------------------------------------------------------------

// returns pow(x,N) only for lanes with nn <= IMAX, and 1.0 otherwise
template <typename T>
template <size_t IMAX> // default should be IMAX=8, but default parameter cannot be used here
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
powLowNLoop(const T xx, const Int_v nn) const
{
  size_t imax = ReduceMax(nn);
  if (imax > IMAX) imax = IMAX;
  Bool_v lowN = (nn <= Int_v(IMAX));

  T res(1.0);
  T nnn(nn);
  for(size_t i = 0; i < imax; ++i) {
    vecCore::MaskedAssign( res, lowN && nnn > 0, res * xx); // only as needed in each lane
    vecCore::MaskedAssign( nnn, lowN && nnn > 0, nnn-1);    // keep track # times / lane
  }
  return res;
}


template <typename T>
VECCORE_ATT_HOST_DEVICE VECCORE_FORCE_INLINE T GXPowVec<T>::
powN(const T x, const Int_v n) const
{
  T xx(x);
  Int_v nn(n);
  vecCore::MaskedAssign( xx, n < 0, T(1.)/x );
  vecCore::MaskedAssign( nn, n < 0, -n );

  // use powLowNLoop() whenever possible
  T result = this->template powLowNLoop<8>(xx,nn);

  Bool_v bigN = (nn > 8);
  if (!vecCore::MaskEmpty(bigN)) {
    const size_t vsize = VectorSize<T>();
    for(size_t i=0; i<vsize; ++i) {
      if( Get(bigN,i) ) Set(result, i, std::pow(xx,nn));
    }
  }
  return result;
}

// -------------------------------------------------------------------

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace
#endif
