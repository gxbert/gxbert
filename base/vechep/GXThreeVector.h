#ifndef GXTHREEVECTOR_H
#define GXTHREEVECTOR_H

/** 
SIMD/SIMT version of HEP3VECTOR of CLHEP
*/

#include "VectorBase.h"
#include "VecHepDefs.h"
//#include "GXGlobal.h"
#include "VecCore/Common.h"

#include <cstdlib>
#include <ostream>
#include <string>

namespace gxbert {
inline namespace GXBERT_IMPL_NAMESPACE {
  //inline namespace VECCORE_IMPL_NAMESPACE {

template <typename T>
class GXThreeVector : VectorBase
{
protected:
  T fx;
  T fy;
  T fz;

public:

  //Constructors
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector() { fx = 0; fy = 0; fz = 0; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector(const T x, const T y, const T z) { fx = x; fy = y; fz = z; }

  //Copy constructor
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector(GXThreeVector<T> const &p) { fx = p.fx; fy = p.fy; fz = p.fz; }

  //Destructor
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  ~GXThreeVector() {}

  //Get methods: x(), y(), z() and GetX(), GetY(), GetZ()
  /*
#define THREEVECTOR_GET(FUNC,VALUE)                \
  VECCORE_ATT_HOST_DEVICE                          \
  VECCORE_FORCE_INLINE                             \          
  T& FUNC () { return VALUE ; }                    \
  VECCORE_ATT_HOST_DEVICE                          \
  VECCORE_FORCE_INLINE                             \                       
  T const& FUNC () const { return VALUE ; } 
  THREEVECTOR_GET(x,fx)
  THREEVECTOR_GET(y,fy)
  THREEVECTOR_GET(z,fz)
  THREEVECTOR_GET(GetX,fx)
  THREEVECTOR_GET(GetY,fy)
  THREEVECTOR_GET(GetX,fz)
#undef THREEVECTOR_GET
  */
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T &x() { return fx; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const &x() const { return fx; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T &y() { return fy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const &y() const { return fy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T &z() { return fz; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const &z() const { return fz; }

  //Set
  VECCORE_ATT_HOST_DEVICE    
  VECCORE_FORCE_INLINE       
  void SetX(T const &x) { fx = x; }

  VECCORE_ATT_HOST_DEVICE    
  VECCORE_FORCE_INLINE       
  void SetY(T const &y) { fy = y; }

 VECCORE_ATT_HOST_DEVICE    
  VECCORE_FORCE_INLINE       
  void SetZ(T const &z) { fz = z; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void Set(T const &x, T const &y, T const &z) { fx = x; fy = y; fz = z; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void Set(const T c) { Set(c, c, c); }

  //Properties
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp2() const { return fx * fx + fy * fy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp() const { return math::Sqrt(Perp2()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Mag2() const { return fx*fx + fy*fy + fz*fz; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Mag() const { return math::Sqrt(Mag2()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Phi() const { return math::ATan2(fy, fx); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Theta() const { return math::ACos(fz / Mag()); }

  //Dot product
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  T Dot(GXThreeVector<T> const &p) const { return fx*p.fx + fy*p.fy + fz*p.fz; }

  //Cross product
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  GXThreeVector<T> Cross(GXThreeVector<T> const &p) const
  {
    return GXThreeVector(fy*p.z()-p.y()*fz, fz*p.x()-p.z()*fx, fx*p.y()-p.x()*fy);
  }

  //Unit : Vector parallel to this, but of length 1.
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXThreeVector<T> Unit() const
  {
    GXThreeVector<T> output(*this);
    const T mag2 = Mag2();
    output /= math::Sqrt(mag2);
    return output;
  }

  // Rotates
  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  GXThreeVector<T>& RotateX(T angle);

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  GXThreeVector<T>& RotateY(T angle);

  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  GXThreeVector<T>& RotateZ(T angle);

  // Rotates reference frame from Uz to newUz (unit vector) (Geant4)
  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  GXThreeVector<T>& RotateUz(GXThreeVector<T>& newUz);

  //operators - index
#define THREEVECTOR_ELEMENT_OP(OPERATOR)                        \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  T operator OPERATOR (int index) const                         \
  {                                                             \
    switch(index) {                                             \
      case 0: return fx;                                        \
      case 1: return fy;                                        \
      case 2: return fz;                                        \
      default: return 0;                                        \
    }                                                           \
  }                                                             \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  T &operator OPERATOR (int index)                              \
  {                                                             \
    switch(index) {                                             \
      case 0: return fx;                                        \
      case 1: return fy;                                        \
      case 2: return fz;                                        \
      default: return 0;                                        \
    }                                                           \
  }
  THREEVECTOR_ELEMENT_OP(())
  THREEVECTOR_ELEMENT_OP([])
#undef THREEVECTOR_ELEMENT_OP

  //operators - assignment
#define THREEVECTOR_ASSIGNMENT_OP(OPERATOR)                     \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  GXThreeVector<T> &operator OPERATOR(const GXThreeVector<T> &p)    \
  {                                                             \
    fx OPERATOR p.fx; 						\
    fy OPERATOR p.fy;                                           \
    fz OPERATOR p.fz;                                           \
    return *this;                                               \
  }                                                              \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  GXThreeVector<T> &operator OPERATOR(const T &c)                 \
  {                                                             \
    fx OPERATOR c; 						\
    fy OPERATOR c;                                              \
    fz OPERATOR c;                                              \
    return *this;                                               \
  }
  THREEVECTOR_ASSIGNMENT_OP(=)
  THREEVECTOR_ASSIGNMENT_OP(+=)
  THREEVECTOR_ASSIGNMENT_OP(-=)
  THREEVECTOR_ASSIGNMENT_OP(*=)
  THREEVECTOR_ASSIGNMENT_OP(/=)
#undef THREEVECTOR_ASSIGNMENT_OP

};

/*
#define THREEVECTOR_ROTATION_AXIS(AXIS, COMPONENT1, COMPONENT2) \
VECCORE_ATT_HOST_DEVICE                                         \
VECCORE_FORCE_INLINE                                            \
template <typename T>                                           \
GXThreeVector<T>& GXThreeVector<T>::Rotate##AXIS(T angle)    	\
{						    		\
  T sinphi = math::Sin(angle);			    		\
  T cosphi = math::Cos(angle);			    		\
  T ta;						    		\
  ta = COMPONENT1 * cosphi - COMPONENT2 * sinphi;  		\
  fz = COMPONENT2 * cosphi + COMPONENT2 * sinphi;               \
  fy = ta;                                                      \
  return *this;                                                 \
}
THREEVECTOR_ROTATION_AXIS(X, fy, fz)
THREEVECTOR_ROTATION_AXIS(Y, fz, fx)
THREEVECTOR_ROTATION_AXIS(Z, fx, fy)
#undef THREEVECTOR_ROTATION_AXIS
*/

template <typename T>
VECCORE_ATT_HOST_DEVICE 
VECCORE_FORCE_INLINE
GXThreeVector<T>& GXThreeVector<T>::RotateUz(GXThreeVector<T>& newUz)
{
  // newUzVector must be normalized !

  T u1 = newUz.x();
  T u2 = newUz.y();
  T u3 = newUz.z();
  T up = math::Sqrt(u1*u1 + u2*u2);

  Mask_v<T> positive = (up > 0.);
  Mask_v<T> negativeZ = (u3 < 0.);

  T invup = Blend(positive,1.0/up, static_cast<T>(0.0));

  T px = fx;
  T py = fy;
  T pz = fz;

  fx = (u1*u3*px - u2*py)*invup + u1*pz;
  fy = (u2*u3*px + u1*py)*invup + u2*pz;
  fz =    -up*px +                u3*pz;

  fx = Blend(positive, fx, Blend(negativeZ, -px, px));
  fy = Blend(positive, fy, Blend(negativeZ,  py, py));
  fz = Blend(positive, fz, Blend(negativeZ, -pz, pz));

  return *this;
}

template <>
VECCORE_ATT_HOST_DEVICE 
VECCORE_FORCE_INLINE
GXThreeVector<double>& GXThreeVector<double>::RotateUz(GXThreeVector<double>& newUz)
{
  //from CLHEP Hep3Vector::rotateUz
  double u1 = newUz.x();
  double u2 = newUz.y();
  double u3 = newUz.z();
  double up = u1*u1 + u2*u2;

  if (up>0) {
    up = std::sqrt(up);
    double px = fx,  py = fy,  pz = fz;
    fx = (u1*u3*px - u2*py)/up + u1*pz;
    fy = (u2*u3*px + u1*py)/up + u2*pz;
    fz =    -up*px +             u3*pz;
  }
  else if (u3 < 0.) { fx = -fx; fz = -fz; }      // phi=0  teta=pi
  else {};
  return *this;
}

#define THREEVECTOR_BINARY_OP(OPERATOR, ASSIGNMENT)                                      \
template <typename T>                                                                    \
VECCORE_FORCE_INLINE                                                                     \
VECCORE_ATT_HOST_DEVICE                                                                  \
GXThreeVector<T> operator OPERATOR(const GXThreeVector<T> &lhs, const GXThreeVector<T> &rhs) \
{                                                                                        \
  GXThreeVector<T> result(lhs);                                                          \
  result ASSIGNMENT rhs;                                                                 \
  return result;                                                                         \
}                                                                                        \
template <typename T, typename ScalarT>                                                  \
VECCORE_FORCE_INLINE                                                                     \
VECCORE_ATT_HOST_DEVICE                                                                  \
GXThreeVector<T> operator OPERATOR(GXThreeVector<T> const &lhs, const ScalarT rhs)       \
{                                                                                        \
  GXThreeVector<T> result(lhs);                                                          \
  result ASSIGNMENT rhs;                                                                 \
  return result;                                                                         \
}                                                                                        \
template <typename T, typename ScalarT>                                                  \
VECCORE_FORCE_INLINE                                                                     \
VECCORE_ATT_HOST_DEVICE                                                                  \
GXThreeVector<T> operator OPERATOR(const ScalarT lhs, GXThreeVector<T> const &rhs)	 \
{                                                                                        \
  GXThreeVector<T> result(lhs);                                                          \
  result ASSIGNMENT rhs;                                                                 \
  return result;                                                                         \
}
THREEVECTOR_BINARY_OP(+, +=)
THREEVECTOR_BINARY_OP(-, -=)
THREEVECTOR_BINARY_OP(*, *=)
THREEVECTOR_BINARY_OP(/, /=)
#undef THREEVECTOR_BINARY_OP

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool operator==(GXThreeVector<T> const &lhs, GXThreeVector<T> const &rhs)
{
  return math::Abs(lhs.x() - rhs.x()) < 0. && 
         math::Abs(lhs.y() - rhs.y()) < 0. && 
         math::Abs(lhs.z() - rhs.z()) < 0. ;
}

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
GXThreeVector<bool> operator!=(GXThreeVector<T> const &lhs, GXThreeVector<T> const &rhs)
{
  return !(lhs == rhs);
}

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
GXThreeVector<T> operator-(GXThreeVector<T> const &v)
{
  return GXThreeVector<T>(-v.x(), -v.y(), -v.z());
}

VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
GXThreeVector<bool> operator!(GXThreeVector<bool> const &v)
{
  return GXThreeVector<bool>(!v.x(), !v.y(), !v.z());
}

template <typename T>
std::ostream &operator<<(std::ostream &os, GXThreeVector<T> const &v)
{
  os << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
  return os;
}

} // end namespace impl
} // end namespace vecCore

#endif
