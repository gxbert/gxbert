#ifndef GXBENCHMARKER_GPU_H
#define GXBENCHMARKER_GPU_H 1

#include "GXTrack.h"
#include "VecHepDefs.h"

namespace gxbert {

  // CUDA

  typedef void (*CudaKernelFunc_t)(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
                                 int nTrackSize, GXTrack *itrack, GXTrack *otrack);

  void CudaTest01(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
                int nTrackSize, GXTrack *itrack, GXTrack *otrack);

  void CudaTest02(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
                int nTrackSize, GXTrack *itrack, GXTrack *otrack);

  CudaKernelFunc_t CudaKernelFunc[] = {
    CudaTest01,
    CudaTest02
  };

} // end namespace gxbert

#endif
