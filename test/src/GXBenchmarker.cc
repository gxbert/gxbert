#include "GXBenchmarker.h"
#include "GXBenchmarker_cpu.h"
#include "GXTrackHandler.h"
#include "VecHepDefs.h"

namespace gxbert {

GXBenchmarker::GXBenchmarker()
    : fNtracks(4992), fRepetitions(1), fVerbosity(1), 
      fMinP(1E+4), fMaxP(1E+4)
{
  fTrackHandler = new GXTrackHandler();
}

GXBenchmarker::~GXBenchmarker()
{
  delete fTrackHandler;
}

int GXBenchmarker::RunBenchmark()
{
  printf(" Run GXBenchmarker with arguments: %d %d %f %f\n", 
    fNtracks, fRepetitions, fMinP, fMaxP);
  printf(" Ntracks       = %d\n", fNtracks);
  printf(" NRepetitions  = %d\n", fRepetitions);
  printf(" MinP (MeV)    = %f\n", fMinP);
  printf(" MaxP (MeV)    = %f\n", fMaxP);

  int errorcode = 0;
  errorcode += RunBenchmarkTests();
  return (errorcode) ? 1 : 0;
}

int GXBenchmarker::RunBenchmarkTests()
{
  //  int mismatches = 0;

  //todo: get number of error from each test and add to mismatches
  RunGeant4();
  RunScalar();
  RunVector();
#ifdef GXBERT_CUDA
  RunCuda();
#endif

//  return mismatches;
  return 0;
}

void GXBenchmarker::RunGeant4()
{
  GXTrack *itrack_aos = (GXTrack *)malloc(fNtracks * sizeof(GXTrack));
  GXTrack *otrack_aos = (GXTrack *)malloc(fNtracks * sizeof(GXTrack));

  double result;
  double elapsedTotal[kNumberTest];
  double elapsedT[kNumberTest];
  double resultTotal[kNumberTest];

  for (int k = 0; k < kNumberTest; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }

  // prepare input tracks
  fTrackHandler->SetRandomStream(1);
  fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
  GXTrack *track_aos = fTrackHandler->GetAoSTracks();

  for (unsigned r = 0; r < fRepetitions; ++r) {
    for (int k = 0; k < kNumberTest; ++k) {
      fTrackHandler->CopyAoSTracks(track_aos, itrack_aos, fNtracks);
      elapsedT[k] = 0.0;
      result = 0;
      elapsedT[k] = Geant4KernelFunc[k](fNtracks, itrack_aos, otrack_aos, result);
      elapsedTotal[k] += elapsedT[k];
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberTest; ++k) {
    printf("%s  Geant4 Total time of %3d reps = %6.3f msec result = %6.3f\n", 
	   TestName[k], fRepetitions, elapsedTotal[k]*1E-6, resultTotal[k]);
  }

  free(itrack_aos);
  free(otrack_aos);
}

void GXBenchmarker::RunScalar()
{
  GXTrack *itrack_aos = (GXTrack *)malloc(fNtracks * sizeof(GXTrack));
  GXTrack *otrack_aos = (GXTrack *)malloc(fNtracks * sizeof(GXTrack));

  double result;
  double elapsedTotal[kNumberTest];
  double elapsedT[kNumberTest];
  double resultTotal[kNumberTest];

  for (int k = 0; k < kNumberTest; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }

    // prepare input tracks
  fTrackHandler->SetRandomStream(1);
  fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
  GXTrack *track_aos = fTrackHandler->GetAoSTracks();

  for (unsigned r = 0; r < fRepetitions; ++r) {
    for (int k = 0; k < kNumberTest; ++k) {
      fTrackHandler->CopyAoSTracks(track_aos, itrack_aos, fNtracks);
      elapsedT[k] = 0.0;
      result = 0;
      elapsedT[k] = ScalarKernelFunc[k](fNtracks, itrack_aos, otrack_aos, result);
      elapsedTotal[k] += elapsedT[k];
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberTest; ++k) {
    printf("%s  Scalar Total time of %3d reps = %6.3f msec result = %6.3f\n",
	   TestName[k], fRepetitions, elapsedTotal[k]*1E-6, resultTotal[k]);
  }

  free(itrack_aos);
  free(otrack_aos);
}

void GXBenchmarker::RunVector()
{
  //  GXTrack *itrack_soa = (GXTrack *)malloc(fNtracks * sizeof(GXTrack));
  //  GXTrack *otrack_soa = (GXTrack *)malloc(fNtracks * sizeof(GXTrack));

  // input SOA tracks
  GXTrackHandler *handler_in = new GXTrackHandler(fNtracks);
  GXTrack_v itrack_soa       = handler_in->GetSoATracks();
  GXTrack_v otrack_soa = handler_in->GetSoATracks();

  // output SOA tracks
  //  GXTrackHandler *handler_out = new GXTrackHandler(fNtracks);

  double result;
  double elapsedTotal[kNumberTest];
  double elapsedT[kNumberTest];
  double resultTotal[kNumberTest];
  for (int k = 0; k < kNumberTest; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }

  // prepare input tracks
  fTrackHandler->SetRandomStream(1);
  fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
  GXTrack_v &track_soa = fTrackHandler->GetSoATracks();

  for (unsigned r = 0; r < fRepetitions; ++r) {
    for (int k = 0; k < kNumberTest; ++k) {
      fTrackHandler->CopySoATracks(track_soa, itrack_soa, fNtracks);
      elapsedT[k] = 0.0;
      result = 0;
      elapsedT[k] = VectorKernelFunc[k](itrack_soa, otrack_soa, result);
      elapsedTotal[k] += elapsedT[k];
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberTest; ++k) {
    printf("%s  Vector Total time of %3d reps = %6.3f msec result = %6.3f\n", 
	   TestName[k], fRepetitions, elapsedTotal[k]*1E-6, resultTotal[k]);
  }

  delete handler_in;
  //  delete handler_out;
}

} // end namespace gxbert
