#include "GXBenchmarker.hh"
#include "GXBenchmarker_gpu.hh"

#include "VecHepDefs.hh"
#include "GXTrackHandler.hh"
#include "GXTrack.hh"
#include "GXCurand.hh"

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
  double resultTotal[kNumberTest];

  double *result_d;

  double *result_h = (double*) calloc(theNBlocks*theNThreads, sizeof(double));
  // malloc(fNtracks*sizeof(double));
  cudaMalloc((void**)&result_d,theNBlocks*theNThreads*sizeof(double));

  for (int k = 0; k < kNumberTest; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }

  fTrackHandler->SetRandomStream(1);
  fTrackHandler->GenerateRandomTracks(fNtracks,fMinP, fMaxP);
  GXTrack* track_aos = fTrackHandler->GetAoSTracks();

  for (unsigned r = 0; r < fRepetitions; ++r) {

    for (int k = 0; k < kNumberTest; ++k) {

      fTrackHandler->CopyAoSTracks(track_aos,itrack_aos,fNtracks);
      cudaMemcpy(itrack_d, track_aos, fNtracks*sizeof(GXTrack), cudaMemcpyHostToDevice);
      
      elapsedEventTime[k] = 0.0;
      
      if(cudaEnabled) {
        cudaEventRecord (start,0);

        //call CUDA kernels
        CudaKernelFunc[k](theNBlocks, theNThreads, randomStates,
                          fNtracks, itrack_d, otrack_d, result_d);

        cudaEventRecord (stop,0);
        cudaEventSynchronize (stop);
        cudaEventElapsedTime (&elapsedEventTime[k],start,stop);

        //copy the result for varification
        cudaMemcpy(result_h,result_d,theNBlocks*theNThreads*sizeof(double),cudaMemcpyDeviceToHost);
        for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) resultTotal[k] += result_h[i];

      }
      elapsedTotal[k] += elapsedEventTime[k]; //ms
      
      cudaMemcpy(itrack_aos, itrack_d, fNtracks*sizeof(GXTrack), cudaMemcpyDeviceToHost);
      cudaMemcpy(otrack_aos, otrack_d, fNtracks*sizeof(GXTrack), cudaMemcpyDeviceToHost);
    }
  }

  for (int k = 0; k < kNumberTest; ++k) {
    printf("%s  Cuda   Total time of %3d reps = %7.4f msec result = %6.3f\n", 
      TestName[k], fRepetitions, elapsedTotal[k], resultTotal[k]);
  }

  //clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(randomStates);
  cudaFree(itrack_d);
  cudaFree(otrack_d);
  cudaFree(result_d);
  free(result_h);

  free(itrack_aos);
  free(otrack_aos);
}

} // end of gxbert namespace
