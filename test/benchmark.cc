#include <iostream>

#include "GXBenchmarker.h"

using namespace gxbert;

int main(int argc, char* argv[])
{
  //default run
  int ntracks = 4992;
  int nrepetitions = 100;
  double minEnergy =  1E+4; // MeV
  double maxEnergy =  minEnergy;

  if(argc >= 2) ntracks =      atoi(argv[1]);
  if(argc >= 3) nrepetitions = atoi(argv[2]);
  if(argc >= 4) {
    minEnergy  =  atof(argv[3]);
    std::cout << "  Min energy = " << minEnergy << std::endl;
    maxEnergy = minEnergy; 
    if(argc >= 5) {
      maxEnergy  =  atof(argv[4]);
    }
    std::cout << "  Max energy = " << maxEnergy << std::endl;
  } 
  else {
    std::cout << "  Min energy = " << minEnergy << std::endl;
    std::cout << "  Max energy = " << maxEnergy << std::endl;
  }

  GXBenchmarker tester;
  tester.SetNTracks(ntracks);
  tester.SetRepetitions(nrepetitions);
  tester.SetMinP(minEnergy);
  tester.SetMaxP(maxEnergy);

  int status = tester.RunBenchmark();

  if(status==1) std::cout << "RunBenchmark Failed" << std::endl;
  return status;
}
