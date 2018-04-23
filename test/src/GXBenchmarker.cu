#include "GXBenchmarker.h"
#include "GXBenchmarker_gpu.h"

#include "VecHepDefs.h"
#include "GXTrackHandler.h"
#include "GXTrack.h"
#include "GXCurand.h"

namespace gxbert {

void GXBenchmarker::RunCuda()
{
  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);
  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  }
  else {
    printf("Waning: No Cuda Capable Device ...");
  }

  //cuda event timing
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  GXTrack* itrack_aos = (GXTrack*) malloc(fNtracks*sizeof(GXTrack));
  GXTrack* otrack_aos = (GXTrack*) malloc(fNtracks*sizeof(GXTrack));

  //allocate memory for input/output tracks
  GXTrack *itrack_d;
  GXTrack *otrack_d;

  cudaMalloc((void**)&itrack_d, fNtracks*sizeof(GXTrack));
  cudaMalloc((void**)&otrack_d, fNtracks*sizeof(GXTrack));

  //set the default number of threads and thread blocks - should be setable
  int theNBlocks  =  26;
  int theNThreads = 192;

  //prepare random engines on the device
  Random_t* randomStates = 0;
  cudaMalloc(&randomStates, theNBlocks*theNThreads* sizeof(curandState));
  GXCurand_Init(randomStates, time(NULL), theNBlocks, theNThreads);

  float elapsedTotal[kNumberTest];
  float elapsedEventTime[kNumberTest];

  for (int k = 0; k < kNumberTest; ++k) elapsedTotal[k] = 0.;

  for (unsigned r = 0; r < fRepetitions; ++r) {

    fTrackHandler->GenerateRandomTracks(fNtracks,fMinP, fMaxP);

    GXTrack* track_aos = fTrackHandler->GetAoSTracks();
    fTrackHandler->SortAoSTracksByEnergy(track_aos,fNtracks);

    for (int k = 0; k < kNumberTest; ++k) {

      fTrackHandler->CopyAoSTracks(track_aos,itrack_aos,fNtracks);
      cudaMemcpy(itrack_d, track_aos, fNtracks*sizeof(GXTrack), cudaMemcpyHostToDevice);
      
      elapsedEventTime[k] = 0.0;
      
      if(cudaEnabled) {
        cudaEventRecord (start,0);

        //call CUDA kernels
        CudaKernelFunc[k](theNBlocks, theNThreads, randomStates,
                          fNtracks, itrack_d, otrack_d);

        cudaEventRecord (stop,0);
        cudaEventSynchronize (stop);
        cudaEventElapsedTime (&elapsedEventTime[k],start,stop);
      }
      elapsedTotal[k] += elapsedEventTime[k]/1000.; //ms to second
      
      cudaMemcpy(itrack_aos, itrack_d, fNtracks*sizeof(GXTrack), cudaMemcpyDeviceToHost);
      cudaMemcpy(otrack_aos, otrack_d, fNtracks*sizeof(GXTrack), cudaMemcpyDeviceToHost);
    }
  }

  for (int k = 0; k < kNumberTest; ++k) {
    printf("%s  Cuda Total time of %3d reps = %7.4f sec\n", 
      TestName[k], fRepetitions, elapsedTotal[k]);
  }

  //clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(randomStates);
  cudaFree(itrack_d);
  cudaFree(otrack_d);

  free(itrack_aos);
  free(otrack_aos);
}

} // end of gxbert namespace
