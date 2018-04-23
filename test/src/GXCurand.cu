#include "GXCurand.h"

#include <iostream>

namespace gxbert {
inline namespace cuda {

__global__
void GXCurand_Init_Kernel(Random_t *randomStates, unsigned long seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &randomStates[tid]);
}

} // end of cuda namespace

bool GXCurand_Init(Random_t *randomStates, unsigned long seed, 
                   int blocksPerGrid, int threadsPerBlock) 
{
  int kstatus = 0;

  gxbert::cuda::GXCurand_Init_Kernel<<<blocksPerGrid,threadsPerBlock>>>
                                                     (randomStates,seed);

  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err ) {
    fprintf(stderr,"GXCurand_Init cudaCheckError() failed at %s : %i : %s\n",
	    __FILE__, __LINE__, cudaGetErrorString(err));
    exit(-1);
  }

  kstatus = cudaThreadSynchronize();
  if (kstatus) std::cout << "GXCurand_Init status = " << kstatus << "\n";
  
  return true;
}

} // end of gxbert namespace
