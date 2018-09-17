#include "timer.h"

#include "GXThreeVector.hh"
#include "LorentzVector.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "GXTrack.hh"
#include "VecHepDefs.hh"

namespace gxbert {

double Geant4Test01(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos,
                    double &result)

{
  CLHEP::Hep3Vector vec3a;
  CLHEP::Hep3Vector vec3b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3b.set(itrack_aos[i].px,itrack_aos[i].py,itrack_aos[i].pz);
    result += vec3a.dot(vec3b);
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double Geant4Test02(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos,
                    double &result)
{
  CLHEP::Hep3Vector vec3a;
  CLHEP::Hep3Vector vec3b;
  CLHEP::Hep3Vector vec3c;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3b.set(itrack_aos[i].px,itrack_aos[i].py,itrack_aos[i].pz);
    vec3c = vec3b.rotateUz(vec3a.unit());
    result += vec3c.perp2();
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double Geant4Boost(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos, 
                   double &result)
{
  CLHEP::HepLorentzVector vec4a;
  CLHEP::HepLorentzVector vec4b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec4a.set(itrack_aos[i].px,itrack_aos[i].py,itrack_aos[i].pz,
              itrack_aos[i].E + itrack_aos[i].m);
    vec4b = vec4a.boost(vec4a.boostVector());
    result += vec4b.perp2();
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

// Scalar

double ScalarTest01(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos,
                    double &result)
{
  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.Set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3b.Set(itrack_aos[i].px,itrack_aos[i].py,itrack_aos[i].pz);
    result += vec3a.Dot(vec3b);
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double ScalarTest02(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos,
                    double &result)
{
  GXThreeVector<double> vec3a;
  GXThreeVector<double> vec3b;
  GXThreeVector<double> vec3c;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec3a.Set(itrack_aos[i].x,itrack_aos[i].y,itrack_aos[i].z);
    vec3b.Set(itrack_aos[i].px,itrack_aos[i].py,itrack_aos[i].pz);
    vec3c = vec3b.RotateUz(vec3a.Unit());
    result += vec3c.Perp2();
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

double ScalarBoost(int ntracks, GXTrack *itrack_aos, GXTrack *otrack_aos,
                   double &result)
{
  LorentzVector<double> vec4a;
  LorentzVector<double> vec4b;
  GXThreeVector<double> vec3;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    vec4a.Set(itrack_aos[i].px,itrack_aos[i].py,itrack_aos[i].pz,
              itrack_aos[i].E + itrack_aos[i].m);
    vec4b = vec4a.Boost(vec4a.BoostVector());
    result +=  vec4b.Perp2();
  }

  elapsedTime = timer.Elapsed();

  return elapsedTime;
}

// Vector

double VectorTest01(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa, double &result)
{
  using Double_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Double_v>();
  int nloop = itrack_soa.size/vsize;

  GXThreeVector<Double_v> vec3a;
  GXThreeVector<Double_v> vec3b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  Double_v result_v(0.);

  timer.Start();

  int ibase = 0;

  Double_v x;
  Double_v y;
  Double_v z;
  Double_v px;
  Double_v py;
  Double_v pz;

  for (int i = 0; i < nloop ; ++i) {
    Load(x, &itrack_soa.x[ibase]);
    Load(y, &itrack_soa.y[ibase]);
    Load(z, &itrack_soa.z[ibase]);
    Load(px, &itrack_soa.px[ibase]);
    Load(py, &itrack_soa.py[ibase]);
    Load(pz, &itrack_soa.pz[ibase]);

    vec3a.Set(x,y,z);
    vec3b.Set(px,py,pz);
    result_v += vec3a.Dot(vec3b);
    ibase += vsize;
  }

  elapsedTime = timer.Elapsed();
  for(int i = 0 ; i < vsize ; i++) result += result_v[i];

  return elapsedTime;
}

double VectorTest02(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa, double &result)
{
  using Double_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Double_v>();
  int nloop = itrack_soa.size/vsize;

  GXThreeVector<Double_v> vec3a;
  GXThreeVector<Double_v> vec3b;
  GXThreeVector<Double_v> vec3c;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  Double_v result_v(0.);
  timer.Start();

  int ibase = 0;

  Double_v x;
  Double_v y;
  Double_v z;
  Double_v px;
  Double_v py;
  Double_v pz;

  for (int i = 0; i < nloop ; ++i) {
    Load(x, &itrack_soa.x[ibase]);
    Load(y, &itrack_soa.y[ibase]);
    Load(z, &itrack_soa.z[ibase]);
    Load(px, &itrack_soa.px[ibase]);
    Load(py, &itrack_soa.py[ibase]);
    Load(pz, &itrack_soa.pz[ibase]);

    vec3a.Set(x,y,z);
    vec3b.Set(px,py,pz);
    vec3c = vec3b.RotateUz(vec3a.Unit());
    result_v += vec3c.Perp2();
    ibase += vsize;
  }
  elapsedTime = timer.Elapsed();
  for(int i = 0 ; i < vsize ; i++) result += result_v[i];

  return elapsedTime;
}

double VectorBoost(GXTrack_v &itrack_soa, GXTrack_v &otrack_soa,
                   double &result)
{
  using Double_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Double_v>();
  int nloop = itrack_soa.size/vsize;

  LorentzVector<Double_v> vec4a;
  LorentzVector<Double_v> vec4b;

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;
  Double_v result_v(0.);

  timer.Start();

  int ibase = 0;

  Double_v px;
  Double_v py;
  Double_v pz;
  Double_v E;
  Double_v m;

  for (int i = 0; i < nloop ; ++i) {
    Load(px, &itrack_soa.px[ibase]);
    Load(py, &itrack_soa.py[ibase]);
    Load(pz, &itrack_soa.pz[ibase]);
    Load(E, &itrack_soa.E[ibase]);
    Load(m, &itrack_soa.m[ibase]);

    vec4a.Set(px, py, pz, E + m);
    vec4b = vec4a.Boost(vec4a.BoostVector());
    result_v += vec4b.Perp2(); 

    //    std::cout << result_v << std::endl;
    ibase += vsize;
  }

  elapsedTime = timer.Elapsed();
  for(int i = 0 ; i < vsize ; i++) result += result_v[i];

  return elapsedTime;
}

} // end namespace gxbert
