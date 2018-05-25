#ifndef GXBENCHMARKER_GPU_H
#define GXBENCHMARKER_GPU_H 1

#include "GXTrack.hh"
#include "VecRng/RngDefs.h"

namespace gxbert {

  // CUDA

  typedef void (*CudaKernelFunc_t)(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
				   int nTrackSize, GXTrack *itrack, GXTrack *otrack, double *result);

  void CudaTest01(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
		  int nTrackSize, GXTrack *itrack, GXTrack *otrack, double *result);

  void CudaTest02(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
		  int nTrackSize, GXTrack *itrack, GXTrack *otrack, double *result);

  void CudaBoost(int blocksPerGrid, int threadsPerBlock, Random_t *devStates,
		 int nTrackSize, GXTrack *itrack, GXTrack *otrack, double *result);

  CudaKernelFunc_t CudaKernelFunc[] = {
    CudaTest01,
    CudaTest02,
    CudaBoost
  };

} // end namespace gxbert

#endif
