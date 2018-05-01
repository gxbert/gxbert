#include "GXThreeVector.h"
#include "LorentzVector.h"
#include "GXTrack.h"

#include "VecRng/RngDefs.h"

namespace gxbert {
inline namespace cuda {

__global__
void KernelTest01(Random_t* devStates, 
                  int nTrackSize, GXTrack* itrack, GXTrack* otrack, double *result)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b;

 __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nTrackSize) {
    vec3a.Set(itrack[tid].x,itrack[tid].y,itrack[tid].z);
    vec3b.Set(itrack[tid].px,itrack[tid].py,itrack[tid].pz);
    tmp += vec3a.Dot(vec3b);
    tid += blockDim.x * gridDim.x;
  }

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];
}

__global__
void KernelTest02(Random_t* devStates, 
                  int nTrackSize, GXTrack* itrack, GXTrack* otrack, double *result)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b;
  GXThreeVector<double> vec3c;

  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nTrackSize) {
    vec3a.Set(itrack[tid].x,itrack[tid].y,itrack[tid].z);
    vec3b.Set(itrack[tid].px,itrack[tid].py,itrack[tid].pz);
    vec3c = vec3b.RotateUz(vec3a.Unit());
    tmp += vec3c.Perp2();
    tid += blockDim.x * gridDim.x;
  }

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];

}

__global__
void KernelBoost(Random_t* devStates, 
                  int nTrackSize, GXTrack* itrack, GXTrack* otrack, double *result)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  LorentzVector<double> vec4a;
  LorentzVector<double> vec4b;
  GXThreeVector<double> vec3;

  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nTrackSize) {
    vec4a.Set(itrack[tid].px,itrack[tid].py,itrack[tid].pz,
              itrack[tid].E + itrack[tid].m);
    vec4b = vec4a.Boost(vec4a.BoostVector());
    tmp +=  vec4b.Perp2();
    tid += blockDim.x * gridDim.x;
  }

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];

}

} // end namespace cuda

// Cuda wrapper

void CudaTest01(int blocksPerGrid, int threadsPerBlock, Random_t* devStates,
 		int nTrackSize, GXTrack* itrack, GXTrack* otrack, double *result) 
{
  gxbert::cuda::KernelTest01<<<blocksPerGrid, threadsPerBlock>>>(devStates,
                nTrackSize,itrack,otrack, result);
}

void CudaTest02(int blocksPerGrid, int threadsPerBlock, Random_t* devStates,
 		int nTrackSize, GXTrack* itrack,GXTrack* otrack, double *result) 
{
  gxbert::cuda::KernelTest02<<<blocksPerGrid, threadsPerBlock>>>(devStates,
                nTrackSize,itrack,otrack,result);
}

void CudaBoost(int blocksPerGrid, int threadsPerBlock, Random_t* devStates,
 		int nTrackSize, GXTrack* itrack,GXTrack* otrack, double *result) 
{
  gxbert::cuda::KernelBoost<<<blocksPerGrid, threadsPerBlock>>>(devStates,
                nTrackSize,itrack,otrack,result);
}

} // end namespace gxbert
