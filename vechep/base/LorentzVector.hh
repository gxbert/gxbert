#ifndef LORENTZVECTOR_H
#define LORENTZVECTOR_H

/*
SIMD/SIMT version of LORENTZVECTOR of CLHEP
*/
#include "VecHepDefs.hh"

#include "VectorBase.hh"
#include "GXThreeVector.hh" //replace by the new one

#include <cstdlib>
#include <ostream>
#include <string>

namespace gxbert {
inline namespace GXBERT_IMPL_NAMESPACE {

template <typename T>
class LorentzVector : VectorBase
{
private:
  GXThreeVector<T> fp;
  T fE;

public:

  //Constructors
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector(const T x, const T y, const T z, const T t)
  {
    fp.Set(x,y,z); fE = t;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector()
  {
    fp.Set(0,0,0); fE = 0;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector(GXThreeVector<T> p, const T t)
  {
    fp = p ;  fE = t;
  }

  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  LorentzVector(LorentzVector<T> const &rhs)
  {
    fp = rhs.fp; fE = rhs.fE;
  }

  /*
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector &operator=(LorentzVector const &rhs)
  {
    fp = rhs.fp; dE = rhs.fE;
    return *this;
  }
  */

  //Get
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& x() { return fp.x(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& x() const { return fp.x(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& y() { return fp.y(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& y() const { return fp.y(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& z() { return fp.z(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& z() const { return fp.z(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& t() { return fE; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& t() const { return fE; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& E() { return fE; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& E() const { return fE; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& e() { return fE; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& e() const { return fE; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& px() { return fp.x(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& px() const { return fp.x(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& py() { return fp.y(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& py() const { return fp.y(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T& pz() { return fp.z(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const& pz() const { return fp.z(); }

  //Set
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetX(T const &x) { fp.SetX(x); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetY(T const &y) { fp.SetY(y); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetZ(T const &z) { fp.SetZ(z); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetT(T const &t) { fE = t; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetE(T const &E) { fE = E; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void Set(T const &x, T const &y, T const &z, T const &t) { fp.Set(x,y,z); fE = t; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void Set(const T c) { Set(c, c, c); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector<T> Vect() { return fp; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  const GXThreeVector<T> Vect() const { return fp; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetVect(const GXThreeVector<T> &p) { fp = p; }

  // Careful: use of Sqrt() makes it slow.
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetVectMag(const GXThreeVector<T> &spatial, T mass)
  {
    SetVect(spatial);
    SetT(math::Sqrt(mass * mass + spatial.Mag2()));
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetVectM(const GXThreeVector<T> &spatial, T mass)
  {
    SetVectMag(spatial, mass);
  }

  // A faster way to set x,y,z, compared to using r,theta,phi arguments, which usually require an Acos(costh) call to get Theta
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetMagCosThPhi(const T mag, const T costh, const T phi, T mass)
  {
    fp.SetMagCosThPhi(mag, costh, phi);
    SetT(math::Sqrt(mass * mass + fp.Mag2()));
  }

  //Properties
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp2() const { return fp.Perp2(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp() const { return fp.Perp(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Mag2() const { return t()*t()- fp.Mag2(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Mag() const { return math::Sqrt(Mag2()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T M2() const { return t()*t()- fp.Mag2(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Phi() const { return math::ATan2(fp.y(), fp.x()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Theta() const { return math::ACos(fp.z() / Mag()) ; }

  // Rotates
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& RotateX(T angle)
  {
    fp.rotateX(angle);
    return *this;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& RotateY(T angle)
  {
    fp.rotateY(angle);
    return *this;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& RotateZ(T angle)
  {
    fp.rotateZ(angle);
    return *this;
  }

  // Rotates reference frame from Uz to newUz (unit vector) (Geant4)
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& RotateUz(LorentzVector<T>& newUz)
  {
    fp.rotateUz(newUz);
    return *this;
  }

  //boost
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector<T> BoostVector() const
  {
    return fp/fE;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& Boost(T bx, T by, T bz);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& Boost(GXThreeVector<T> b)
  {
    return Boost(b.x(),b.y(),b.z());
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& BoostX(T bbeta);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& BoostY(T bbeta);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector<T>& BoostZ(T bbeta);

  //operators - index
#define LORENTZVECTOR_INDEX_OP(OPERATOR)                        \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  T operator OPERATOR(int index) const                          \
  {                                                             \
    switch(index) {                                             \
      case 0: return fp.fx;                                     \
      case 1: return fp.fy;                                     \
      case 2: return fp.fz;                                     \
      case 3: return fE;                                        \
      default: return 0;                                        \
    }                                                           \
  }                                                             \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  T &operator OPERATOR(int index)                               \
  {                                                             \
    switch(index) {                                             \
      case 0: return fp.fx;                                     \
      case 1: return fp.fy;                                     \
      case 2: return fp.fz;                                     \
      case 3: return fE;                                        \
      default: return 0;                                        \
    }                                                           \
  }
  LORENTZVECTOR_INDEX_OP(())
  LORENTZVECTOR_INDEX_OP([])
#undef LORENTZVECTOR_INDEX_OP

#define LORENTZVECTOR_ASSIGNMENT_OP(OPERATOR)                     \
  VECCORE_ATT_HOST_DEVICE                                         \
  VECCORE_FORCE_INLINE                                            \
  LorentzVector<T> &operator OPERATOR(const LorentzVector<T> &p)  \
  {                                                               \
    fp OPERATOR p.fp;                                             \
    fE OPERATOR p.fE;                                             \
    return *this;                                                 \
  }                                                               \
  template <typename T1>                                          \
  VECCORE_ATT_HOST_DEVICE                                         \
  VECCORE_FORCE_INLINE                                            \
  LorentzVector<T> &operator OPERATOR(const LorentzVector<T1> &p) \
  {                                                               \
    fp OPERATOR p.fp;                                             \
    fE OPERATOR p.fE;                                             \
    return *this;                                                 \
  }                                                               \
  VECCORE_ATT_HOST_DEVICE                                         \
  VECCORE_FORCE_INLINE                                            \
  LorentzVector<T> &operator OPERATOR(const T &c)                 \
  {                                                               \
    fp OPERATOR c;                                                \
    fE OPERATOR c;                                                \
    return *this;                                                 \
  }
  LORENTZVECTOR_ASSIGNMENT_OP(=)
  LORENTZVECTOR_ASSIGNMENT_OP(+=)
  LORENTZVECTOR_ASSIGNMENT_OP(-=)
  LORENTZVECTOR_ASSIGNMENT_OP(*=)
  LORENTZVECTOR_ASSIGNMENT_OP(/=)
#undef LORENTZVECTOR_ASSIGNMENT_OP
};

template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
LorentzVector<T>& LorentzVector<T>::Boost(T bx, T by, T bz)
{
  const T kZero(0.0);
  const T kOne(1.0);
  const T b2 = bx*bx + by*by + bz*bz;
  const T ggamma = kOne / math::Sqrt(kOne - b2);
  const T bp = bx*x() + by*y() + bz*z();

  Mask_v<T> positive = (b2 > kZero);
  const T gamma2 =  Blend(positive,(ggamma - kOne)/b2, kZero);

  SetX(x() + gamma2*bp*bx + ggamma*bx*t());
  SetY(y() + gamma2*bp*by + ggamma*by*t());
  SetZ(z() + gamma2*bp*bz + ggamma*bz*t());
  SetT(ggamma*(t() + bp));

  return *this;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
LorentzVector<T>& LorentzVector<T>::BoostX(T bbeta)
{
  T b2 = bbeta*bbeta;

  //check beta >= 1 (speed of light): no boost done along x
  Mask_v<T> ge_c = (b2 >= 1);
  bbeta = Blend(ge_c, 0.0, bbeta);
  T ggamma = Blend(ge_c, 1.0, math::Sqrt(1./(1-b2)));

  T tt = fE;
  fE = ggamma*(fE + bbeta*fp.getX());
  fp.SetX(ggamma*(fp.getX() + bbeta*tt));

  return *this;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
LorentzVector<T>& LorentzVector<T>::BoostY(T bbeta)
{
  T b2 = bbeta*bbeta;

  //check beta >= 1 (speed of light): no boost done along y
  Mask_v<T> ge_c = (b2 >= 1);
  bbeta = Blend(ge_c, 0.0, bbeta);
  T ggamma = Blend(ge_c, 1.0, math::Sqrt(1./(1-b2)));

  T tt = fE;
  fE = ggamma*(fE + bbeta*fp.GetY());
  fp.SetY(ggamma*(fp.GetY() + bbeta*tt));

  return *this;
}

template <typename T>
VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
LorentzVector<T>& LorentzVector<T>::BoostZ(T bbeta)
{
  const T kZero(0.0);
  const T kOne(1.0);
  const T b2 = bbeta*bbeta;

  //check beta >= 1 (speed of light): no boost done along z
  Mask_v<T> ge_c = (b2 >= 1);
  bbeta = Blend(ge_c, kZero, bbeta);
  const T ggamma = Blend(ge_c, kOne, math::Sqrt(kOne / (kOne - b2)));

  const T tt = fE;
  fE = ggamma * (fE + bbeta*fp.GetZ());
  fp.SetZ(ggamma*(fp.GetZ() + bbeta*tt));

  return *this;
}


#define LORENTZVECTOR_BINARY_OP(OPERATOR, ASSIGNMENT)                                         \
template <typename T, typename OtherT>                                                        \
VECCORE_FORCE_INLINE                                                                          \
VECCORE_ATT_HOST_DEVICE                                                                       \
LorentzVector<T> operator OPERATOR(const LorentzVector<T> &lhs, const LorentzVector<OtherT> &rhs)  \
{                                                                                             \
  LorentzVector<T> result(lhs);                                                               \
  result ASSIGNMENT rhs;                                                                      \
  return result;                                                                              \
}                                                                                             \
template <typename T, typename ScalarT>                                                       \
VECCORE_FORCE_INLINE                                                                          \
VECCORE_ATT_HOST_DEVICE                                                                       \
LorentzVector<T> operator OPERATOR(LorentzVector<T> const &lhs, const ScalarT rhs)            \
{                                                                                             \
  LorentzVector<T> result(lhs);                                                               \
  result ASSIGNMENT rhs;                                                                      \
  return result;                                                                              \
}                                                                                             \
template <typename T, typename ScalarT>                                                       \
VECCORE_FORCE_INLINE                                                                          \
VECCORE_ATT_HOST_DEVICE                                                                       \
LorentzVector<T> operator OPERATOR(const ScalarT lhs, LorentzVector<T> const &rhs)	      \
{                                                                                             \
  LorentzVector<T> result(lhs);                                                               \
  result ASSIGNMENT rhs;                                                                      \
  return result;                                                                              \
}
LORENTZVECTOR_BINARY_OP(+, +=)
LORENTZVECTOR_BINARY_OP(-, -=)
LORENTZVECTOR_BINARY_OP(*, *=)
LORENTZVECTOR_BINARY_OP(/, /=)
#undef LORENTZVECTOR_BINARY_OP

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool operator==(LorentzVector<T> const &lhs, LorentzVector<T> const &rhs)
{
  T kZero(0.0);
  return math::Abs(lhs.x() - rhs.x()) < kZero &&
         math::Abs(lhs.y() - rhs.y()) < kZero &&
         math::Abs(lhs.z() - rhs.z()) < kZero &&
         math::Abs(lhs.t() - rhs.t()) < kZero ;
}

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
LorentzVector<bool> operator!=(LorentzVector<T> const &lhs, LorentzVector<T> const &rhs)
{
  return !(lhs == rhs);
}

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
LorentzVector<T> operator-(LorentzVector<T> const &v)
{
  return LorentzVector<T>(-v.x(), -v.y(), -v.z(), -v.t());
}

VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
LorentzVector<bool> operator!(LorentzVector<bool> const &v)
{
  return LorentzVector<bool>(!v.x(), !v.y(), !v.z(), !v.t());
}

template <typename T>
std::ostream &operator<<(std::ostream &os, LorentzVector<T> const &v)
{
  os << "(" << v.x() << ", " << v.y() << ", "  << v.z() << ", " << v.t() << ")";
  return os;
}

} // end namespace impl
} // end namespace vecCore

#endif
