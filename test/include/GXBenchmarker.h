#ifndef GXBenchmarker_H
#define GXBenchmarker_H 1

#include "VecHepDefs.h"

namespace gxbert {

enum TestIndex { 
  kNullTest = -1, 
  kTest01, 
  kTest02, 
  kTestBoost, 
  kNumberTest 
};

static const char *TestName[kNumberTest] = {
  "ThreeVectorDot     ", 
  "ThreeVectorRotateUz",
  "LoentzVectorBoost  "
};


class GXTrackHandler;

class GXBenchmarker {

public:
  GXBenchmarker();
  ~GXBenchmarker();

  int RunBenchmark();

  void SetNTracks(const int ntracks) { fNtracks = ntracks; }
  void SetRepetitions(const unsigned repetitions) { fRepetitions = repetitions; }

  void SetMinP(double pMin) { fMinP = pMin; }
  void SetMaxP(double pMax) { fMaxP = pMax; }

private:
  int RunBenchmarkTests();

  void RunGeant4();
  void RunScalar();
  void RunVector();

#ifdef GXBERT_CUDA
  void RunCuda();
#endif

private:

  GXTrackHandler *fTrackHandler;

  int fNtracks;
  unsigned fRepetitions;
  int fVerbosity;

  double fMinP;
  double fMaxP;
};

} // end namespace gxbert

#endif
