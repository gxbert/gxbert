// @File: ApproxEqual.hh
//
// @Purpose: Convenience functions for FP comparisons
//
// 20180601 Guilherme Lima -- Vectorized version, extended original from vecgeom/test/unit_test/ApproxEqual.h


#ifndef APPROXEQUAL_HH
#define APPROXEQUAL_HH

#include <cmath>
#include "VecCore/VecCore"

constexpr double kApproxEqualTolerance = 1.E-6;
constexpr double kInfinity = vecCore::NumericLimits<double>::Max();

// Return true if the double x is approximately equal to y
//
// Process:
//
// Return true is x if less than kApproxEqualTolerance from y

template <typename T>
vecCore::Mask_v<T> ApproxEqual(T const& x, T const& y)
{
  using vecCore::math::Abs;
  vecCore::Mask_v<T> result( x == y );

  if ( vecCore::VectorSize<T>() == 1 ) {
    if ( vecCore::MaskFull( result ) ) return result;
  }

  //  } else if (x * y == 0.0)
  //    T diff = Abs(x - y);
  //    return diff < kApproxEqualTolerance;
  //  }
  //result = result | 
  //else if (Abs(x) > kInfinity || Abs(y) > kInfinity) {
  //  // handle comparisons to infinity
  //  return (x * y > 0) && Abs(x) > kInfinity && Abs(y) > kInfinity;
  //} else {

  T diff  = Abs(x - y);
  //T abs_x = Abs(x), abs_y = Abs(y);
  result = result | (diff < kApproxEqualTolerance);

  //.. handle comparisons to infinity
  //return (x * y > 0) && Abs(x) > kInfinity && Abs(y) > kInfinity;
  result = result | ( (x * y > T(0.)) & (Abs(x) > kInfinity) & (Abs(y) > kInfinity) );

  return result;
}

// Return true if the 3vector check is approximately equal to target
template <typename Vec_t>
bool ApproxEqualVec(const Vec_t &check, const Vec_t &target)
{
  return (ApproxEqual(check.x(), target.x()) && ApproxEqual(check.y(), target.y()) &&
          ApproxEqual(check.z(), target.z()))
             ? true
             : false;
}

#endif
