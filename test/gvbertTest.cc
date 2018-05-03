#include <Vc/Vc>
#include "VecCore/Timer.h"
#include "VecRng/RngDefs.h"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"

using namespace gxbert;

int main(int argc, char* argv[]){

  //default run with a fixed target and energy
  double energy =  1.3E+3; // 1.3 GeV
  int event_per_job = 1024*64;

  if(argc >= 2) energy = atof(argv[1]);
  if(argc >= 3) event_per_job = atoi(argv[2]);

  //timer
  Timer<nanoseconds> timer;
  double timeElapsed = 0.;

  //vector size 
  using Double_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Double_v>();

  //prepare input
  GXTrackHandler *handler = new GXTrackHandler(event_per_job);
  handler->GenerateRandomTracks(event_per_job, energy, energy);
  GXTrack_v track_soa = handler->GetSoATracks();

  printf("Number of tracks = %d\n",track_soa.size);

  timer.Start();

  //event loop
  int ibase = 0;

  for (int i = 0; i < event_per_job/vsize ; ++i) {
    //make a reaction - add the top level <VectorBackend> interface here
    // result = cascadeInterface->ApplyYourself<Vector>(projectile, nucleus);

    ibase +=  vsize;
  }

  timeElapsed = timer.Elapsed();
  printf("Time Elapsed = %7.4f ms\n",timeElapsed*1E-6);

  //clean up
  delete handler;

  return 0;
}
