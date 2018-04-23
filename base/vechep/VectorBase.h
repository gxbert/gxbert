#ifndef VECMATH_VECTORBASE_H
#define VECMATH_VECTORBASE_H

#ifndef VECCORE_NVCC
#include <stdlib.h>
#include "VecCore/Common.h"
#endif

namespace vecCore {
inline namespace VECCORE_IMPL_NAMESPACE {

#ifndef VECCORE_NVCC

class VectorBase {
public:

  /*
  VECCORE_FORCE_INLINE
  void *operator new(size_t size) {
    return aligned_malloc(VECCORE_SIMD_ALIGN, size);
  }

  VECCORE_FORCE_INLINE
  void *operator new[](size_t size) {
    return aligned_malloc(VECCORE_SIMD_ALIGN, size);
  }
  */
};
#else
class VectorBase {};
#endif

} // end namespace impl
} // end namespace vecCore

#endif // VECMATH_VECTORBASE_H
