#ifndef GXBENCHMARKER_CPU_H
#define GXBENCHMARKER_CPU_H 1

#include "GXTrack.h"

namespace gxbert {

  // Geant4

  typedef double (*Geant4KernelFunc_t)(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos);

  double Geant4Test01(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos);
  double Geant4Test02(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos);

  Geant4KernelFunc_t Geant4KernelFunc[] = {
    Geant4Test01,  
    Geant4Test02
  };

  // Scalar

  typedef double (*ScalarKernelFunc_t)(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos);

  double ScalarTest01(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos);
  double ScalarTest02(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos);

  ScalarKernelFunc_t ScalarKernelFunc[] = {
    ScalarTest01,  
    ScalarTest02
  };

  // Vector

  typedef double (*VectorKernelFunc_t)(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa);

  double VectorTest01(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa);
  double VectorTest02(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa);

  VectorKernelFunc_t VectorKernelFunc[] = {
    VectorTest01,  
    VectorTest02
  };

} // end namespace gxbert

#endif
