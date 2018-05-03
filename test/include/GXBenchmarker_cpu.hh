#ifndef GXBENCHMARKER_CPU_H
#define GXBENCHMARKER_CPU_H 1

#include "GXTrack.hh"

namespace gxbert {

  // Geant4

  typedef double (*Geant4KernelFunc_t)(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, 
                                       double &result);

  double Geant4Test01(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, double &result);
  double Geant4Test02(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, double &result);
  double Geant4Boost(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, double &result);

  Geant4KernelFunc_t Geant4KernelFunc[] = {
    Geant4Test01,  
    Geant4Test02,  
    Geant4Boost
  };

  // Scalar

  typedef double (*ScalarKernelFunc_t)(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos,
                                        double &result);

  double ScalarTest01(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, double &result);
  double ScalarTest02(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, double &result);
  double ScalarBoost(int ntrack, GXTrack *itrack_aos, GXTrack *otrack_aos, double &result);

  ScalarKernelFunc_t ScalarKernelFunc[] = {
    ScalarTest01,  
    ScalarTest02,
    ScalarBoost
  };

  // Vector

  typedef double (*VectorKernelFunc_t)(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa,
                                        double &result);

  double VectorTest01(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa, double &result);
  double VectorTest02(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa, double &result);
  double VectorBoost(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa, double &result);

  VectorKernelFunc_t VectorKernelFunc[] = {
    VectorTest01,  
    VectorTest02,
    VectorBoost
  };

} // end namespace gxbert

#endif
