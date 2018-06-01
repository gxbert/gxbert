// @File: ApproxEqual.hh
//
// @Purpose: Convenience functions for FP comparisons
//
// 20180601 Guilherme Lima -- Shamelessly copied from VecGeom unit tests

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
bool ApproxEqual(T const& x, T const& y)
{
  using vecCore::math::Abs;
  if (x == y) {
    return true;
  } else if (x * y == 0.0) {
    T diff = Abs(x - y);
    return diff < kApproxEqualTolerance;
  } else if (Abs(x) > kInfinity || Abs(y) > kInfinity) {
    // handle comparisons to infinity
    return (x * y > 0) && Abs(x) > kInfinity && Abs(y) > kInfinity;
  } else {
    double diff  = Abs(x - y);
    double abs_x = Abs(x), abs_y = Abs(y);
    return diff / (abs_x + abs_y) < kApproxEqualTolerance;
  }
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
