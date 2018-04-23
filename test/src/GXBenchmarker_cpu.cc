#include "VecCore/Timer.h"

#include "GXThreeVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "GXTrack.h"
#include "VecHepDefs.h"

namespace gxbert {

double Geant4Test01(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos)
{
  CLHEP::Hep3Vector vec3a;
  CLHEP::Hep3Vector vec3b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3a.dot(vec3b);
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double Geant4Test02(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos)
{
  CLHEP::Hep3Vector vec3a;
  CLHEP::Hep3Vector vec3b(1.,0.,0.);
  CLHEP::Hep3Vector vec3c;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3c = vec3a.rotateUz(vec3b);
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

// Scalar

double ScalarTest01(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos)
{
  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.Set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3a.Dot(vec3b);
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double ScalarTest02(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos)
{
  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b(1.,0.,0.);
  GXThreeVector<double> vec3c;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.Set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3c = vec3a.RotateUz(vec3b);
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

// Vector

double VectorTest01(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa)
{
  using Double_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Double_v>();
  int nloop = itrack_soa.size/vsize;

  GXThreeVector<Double_v> vec3a;
  GXThreeVector<Double_v> vec3b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  int ibase = 0;

  for (int i = 0; i < nloop ; ++i) {
    vec3a.Set(itrack_soa.x[ibase],itrack_soa.y[ibase],itrack_soa.z[ibase]);
    vec3a.Dot(vec3b);
    ibase += vsize;
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double VectorTest02(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa)
{
  using Double_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Double_v>();
  int nloop = itrack_soa.size/vsize;

  GXThreeVector<Double_v> vec3a;
  GXThreeVector<Double_v> vec3b(1.,0.,0.);
  GXThreeVector<Double_v> vec3c;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  int ibase = 0;

  for (int i = 0; i < nloop ; ++i) {
    vec3a.Set(itrack_soa.x[ibase],itrack_soa.y[ibase],itrack_soa.z[ibase]);
    vec3c = vec3a.RotateUz(vec3b);
    ibase += vsize;
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

} // end namespace gxbert
