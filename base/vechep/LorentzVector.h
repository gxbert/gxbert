#ifndef LORENTZVECTOR_H
#define LORENTZVECTOR_H

/* 
SIMD/SIMT version of LORENTZVECTOR of CLHEP
*/

#include "VectorBase.h"
#include "GXThreeVector.h" //replace by the new one

#include <cstdlib>
#include <ostream>
#include <string>

namespace vecCore {
inline namespace VECCORE_IMPL_NAMESPACE {

template <typename T>
class LorentzVector : VectorBase
{
private:
  ThreeVector<T> fp;
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
    fp(0); fE = 0;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector(ThreeVector<T> p, const T t)
  {
    fp = p ;  fE = t;
  }

  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  LorentzVector(LorentzVector<T> const &rhs)
  {
    fp = rhs.fp; fE = rhs.fE;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  LorentzVector &operator=(LorentzVector const &rhs)
  {
    fp = rhs.fp; dE = rhs.fE;
    return *this;
  }

 //Get
  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE                       
  T& x() { return fp.fx; }                                          

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE                       
  T const& x() const { return fp.fx; }                                          

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE                       
  T& y() { return fp.fy; }                                          

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE                       
  T const& y() const { return fp.fy; }                                          

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE                       
  T& z() { return fp.fz; }                                          

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE                       
  T const& z() const { return fp.fz; }                                          

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

  //Set
  VECCORE_ATT_HOST_DEVICE    
  VECCORE_FORCE_INLINE       
  void SetX(T const &x) { fp.fx = x; }

  VECCORE_ATT_HOST_DEVICE    
  VECCORE_FORCE_INLINE       
  void SetY(T const &y) { fp.fy = y; }

  VECCORE_ATT_HOST_DEVICE    
  VECCORE_FORCE_INLINE       
  void SetZ(T const &z) { fp.fz = z; }

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
  void SetVect(const ThreeVector<T> &p) { fp = p; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetVectMag(const ThreeVector<T> &spatial, T magnitude) 
  { 
    SetVect(spatial);
    SetT(math::Sqrt(magnitude * magnitude + spatial * spatial));
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetVectM(const ThreeVector<T> &spatial, T magnitude) 
  { 
    SetVecMag(spatial,magnitude);
  }

  //Properties
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp2() const { return fp.Perp2(); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp() const { return fp.Perp()); }

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
  T Phi() const { return math::ATan2(fp.fy, fp.fx); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Theta() const { return math::ACos(fp.fz / Mag()); }

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
  LorentzVector<T>& Boost(T bx, T by, T bz); 

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  LorentzVector<T>& Boost(const ThreeVector<T>& b) 
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
      case 0: return fpfx;                                      \
      case 1: return fpfy;                                      \
      case 2: return fpfz;                                      \
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
  T b2 = bx*bx + by*by + bz*bz;
  T ggamma = 1.0 / math::Sqrt(1.0 - b2);
  T bp = bx*x() + by*y() + bz*z();

  Mask_v<T> positive = (b2 > 0.);
  T gamma2 =  Blend(positive,(ggamma - 1.0)/b2, static_cast<T>(0.0));

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
  T ggamma = Blend(ge_c, 1.0, math::Sqrt(1./(1-b2));

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
  T ggamma = Blend(ge_c, 1.0, math::Sqrt(1./(1-b2));

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
  T b2 = bbeta*bbeta;

  //check beta >= 1 (speed of light): no boost done along z
  Mask_v<T> ge_c = (b2 >= 1);
  bbeta = Blend(ge_c, 0.0, bbeta);
  T ggamma = Blend(ge_c, 1.0, math::Sqrt(1./(1-b2));

  T tt = fE;
  fE = ggamma*(fE + bbeta*fp.GetZ());
  fp.SetZ(ggamma*(fp.GetZ() + bbeta*tt));

  return *this;
} 


#define LORENTZVECTOR_BINARY_OP(OPERATOR, ASSIGNMENT)                                         \
VECCORE_FORCE_INLINE                                                                          \
VECCORE_ATT_HOST_DEVICE                                                                       \
LorentzVector<T> operator OPERATOR(const LorentzVector<T> &lhs, const LorentzVector<T> &rhs)  \
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

VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool operator==(LorentzVector<Real_t> const &lhs, LorentzVector<Real_t> const &rhs)
{
  return math::Abs(lhs.x() - rhs.x()) < 0. && 
         math::Abs(lhs.y() - rhs.y()) < 0. && 
         math::Abs(lhs.z() - rhs.z()) < 0. && 
         math::Abs(lhs.t() - rhs.t()) < 0. ;
}

VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
LorentzVector<bool> operator!=(LorentzVector<Real_t> const &lhs, LorentzVector<Real_t> const &rhs)
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
