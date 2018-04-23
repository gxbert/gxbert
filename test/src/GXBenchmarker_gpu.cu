#include "GXThreeVector.h"
#include "VecHepDefs.h"
#include "GXTrack.h"

namespace gxbert {
inline namespace cuda {

__global__
void KernelTest01(Random_t* devStates, 
                  int nTrackSize, GXTrack* itrack, GXTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b;

  double result = 0;

  while (tid < nTrackSize) {
    vec3a.Set(itrack[tid].x,itrack[tid].y,itrack[tid].z);
    result += vec3a.Dot(vec3b);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelTest02(Random_t* devStates, 
                  int nTrackSize, GXTrack* itrack, GXTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b(1.,0.,0.);
  GXThreeVector<double> vec3c;

  while (tid < nTrackSize) {
    vec3a.Set(itrack[tid].x,itrack[tid].y,itrack[tid].z);
    vec3c = vec3a.RotateUz(vec3b);
    tid += blockDim.x * gridDim.x;
  }
}

} // end namespace cuda

// Cuda wrapper

void CudaTest01(int blocksPerGrid, int threadsPerBlock, Random_t* devStates,
 		int nTrackSize, GXTrack* itrack, GXTrack* otrack) 
{
  gxbert::cuda::KernelTest01<<<blocksPerGrid, threadsPerBlock>>>(devStates,
                nTrackSize,itrack,otrack);
}

void CudaTest02(int blocksPerGrid, int threadsPerBlock, Random_t* devStates,
 		int nTrackSize, GXTrack* itrack,GXTrack* otrack) 
{
  gxbert::cuda::KernelTest02<<<blocksPerGrid, threadsPerBlock>>>(devStates,
                nTrackSize,itrack,otrack);
}

} // end namespace gxbert
