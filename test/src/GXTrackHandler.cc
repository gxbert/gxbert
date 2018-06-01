#include "GXTrackHandler.hh"

#include "mm_malloc.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

//#include "VecRng/MRG32k3a.h"

namespace gxbert {


GXTrackHandler::GXTrackHandler() 
  : fNumberOfTracks(0), fTrack_aos(0), fBuffer(0), fMass(938.272013) // proton mass by default
{
  // Scalar MRG32k3a
  fRNG = new  vecRng::cxx::MRG32k3a<ScalarBackend>;
  fRNG->Initialize();
}

GXTrackHandler::GXTrackHandler(size_t nTracks)
  : GXTrackHandler()
{
  Allocate(nTracks);
}

GXTrackHandler::~GXTrackHandler() { Deallocate(); }

void GXTrackHandler::SetNumberOfTracks(size_t nTracks) 
{ 
  fNumberOfTracks = nTracks; 
}

void GXTrackHandler::Deallocate()
{
  if (fNumberOfTracks > 0) {
    _mm_free(fTrack_aos);
    _mm_free(fBuffer);
  }
}

void GXTrackHandler::Allocate(size_t nTracks)
{
  const int blockSize = 64; // Bytes

  SetNumberOfTracks(nTracks);

  if (fNumberOfTracks > 0) {

    unsigned int numAlloc = ((fNumberOfTracks - 1) / blockSize + 1) * blockSize;

    unsigned int memSizeAlloc = numAlloc * sizeof(GXTrack);

    // allocation for aos
    fTrack_aos = (GXTrack *)_mm_malloc(memSizeAlloc, blockSize);

    memSizeAlloc += blockSize;

    // allocation for soa
    char *soaBuffer = (char *)_mm_malloc(memSizeAlloc, blockSize);
    fBuffer = soaBuffer;

    // stride for GXTrack_v.numTracks
    soaBuffer += blockSize;   

    //    const int offset_int = numAlloc * sizeof(int);
    const int offset_double = numAlloc * sizeof(double);

    // set ptr to each element of GXTrack_v
    fTrack_soa.x = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.y = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.z = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.t = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.px = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.py = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.pz = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.E = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.q = (double *)soaBuffer;
    soaBuffer += offset_double;

    fTrack_soa.m = (double *)soaBuffer;

    assert((long)soaBuffer + offset_double <= (long)fBuffer + memSizeAlloc);
  }
}

void GXTrackHandler::Reallocate(size_t nTracks)
{
  Deallocate();
  Allocate(nTracks);
}

double GXTrackHandler::Random() 
{
  return fRNG->Uniform<ScalarBackend>();
}

void GXTrackHandler::SetRandomStream(long streamId) 
{
  fRNG->Initialize(streamId);
}

void GXTrackHandler::GenerateRandomTracks(size_t nTracks, double minP, double maxP)
{
  Reallocate(nTracks);

  // constants - move to a header file
  double const pi = 3.14159265358979323846;
  double const ecalRmim = 1290.;
  double const ecalRmax = 1520.;
  double const ecalZmax = 3000.;

  fTrack_soa.size = fNumberOfTracks;

  for (size_t i = 0; i < fNumberOfTracks; ++i) {
    double rho, z, p;
    double theta, phi;
    double cosphi, sinphi;
    double sintheta, tantheta, costheta;

    rho = ecalRmim + (ecalRmax - ecalRmim) * Random();

    if (minP == maxP) {
      p = maxP;
    }
    else {
      do {
        p = minP - 0.2 * (maxP - minP) * std::log(Random());
      } while (p > maxP);
    }

    z = ecalZmax * (2 * Random() - 1.0);
    phi = 2 * pi * Random();
    tantheta = rho / z;
    theta = std::atan(tantheta); // (rho/z);

    cosphi = std::cos(phi);
    sinphi = std::sin(phi);

    costheta = std::cos(theta);
    sintheta = std::sin(theta);

    (fTrack_soa.x)[i] = fTrack_aos[i].x = rho * cosphi;
    (fTrack_soa.y)[i] = fTrack_aos[i].y = rho * sinphi;
    (fTrack_soa.z)[i] = fTrack_aos[i].z = z;
    (fTrack_soa.t)[i] = fTrack_aos[i].t = 0.;

    (fTrack_soa.px)[i] = fTrack_aos[i].px = p * sintheta * cosphi;
    (fTrack_soa.py)[i] = fTrack_aos[i].py = p * sintheta * sinphi;
    (fTrack_soa.pz)[i] = fTrack_aos[i].pz = p * costheta;
    (fTrack_soa.E)[i] = fTrack_aos[i].E = p*p/(sqrt(p*p + fMass*fMass) + fMass);

    (fTrack_soa.q)[i] = fTrack_aos[i].q = 1.0;
    (fTrack_soa.m)[i] = fTrack_aos[i].m = fMass;
  }
}

  // note: momentum range in MeV units!!!
void GXTrackHandler::GenerateTracksAlongSameDirection(size_t nTracks, double *posdir, double minP, double maxP)
{
  // make sure input direction vector is a unit vector
  assert(posdir[3]*posdir[3] + posdir[4]*posdir[4] + posdir[5]*posdir[5] - 1.0 < 1.0e-9);

  Reallocate(nTracks);

  fTrack_soa.size = fNumberOfTracks;
  for (size_t i = 0; i < fNumberOfTracks; ++i) {
    double p = maxP;
    if (minP < maxP) p = minP + (maxP - minP) * std::log(Random());

    (fTrack_soa.x)[i] = fTrack_aos[i].x = posdir[0];
    (fTrack_soa.y)[i] = fTrack_aos[i].y = posdir[1];
    (fTrack_soa.z)[i] = fTrack_aos[i].z = posdir[2];
    (fTrack_soa.t)[i] = fTrack_aos[i].t = 0.;

    (fTrack_soa.px)[i] = fTrack_aos[i].px = p * posdir[3];
    (fTrack_soa.py)[i] = fTrack_aos[i].py = p * posdir[4];
    (fTrack_soa.pz)[i] = fTrack_aos[i].pz = p * posdir[5];
    (fTrack_soa.E)[i] = fTrack_aos[i].E = p*p/(sqrt(p*p + fMass*fMass) + fMass);

    (fTrack_soa.q)[i] = fTrack_aos[i].q = 1.0;
    (fTrack_soa.m)[i] = fTrack_aos[i].m = fMass;
    if (i == 0) {
      double Ekin = fTrack_aos[i].E;
      double Etot = Ekin + fMass;
      std::cerr<<"GXTracksHandler: i="<< i <<'/'<< fNumberOfTracks <<": "<< fTrack_aos[i].x <<" "<< fTrack_aos[i].y <<" "<< fTrack_aos[i].z <<" "<< fTrack_aos[i].px <<" "<< fTrack_aos[i].py <<" "<< fTrack_aos[i].pz  <<" Ek="<< Ekin <<" Etot="<< Etot <<"\n";
    }
  }
}

// utility functions - can be elsewhere

void GXTrackHandler::SortAoSTracksByEnergy(GXTrack *AoS, size_t ntracks)
{
  // sort AoS tracks by energy in ascending order.
  std::sort(AoS, AoS + ntracks, 
    [](GXTrack const &a, GXTrack const &b) { return a.E < b.E; });
}

void GXTrackHandler::SortSoATracksByEnergy(GXTrack_v &SoA, size_t ntracks)
{
  // temporary AoS tracks
  const int blockSize = 64; // Bytes
  GXTrack *AoS = (GXTrack *)_mm_malloc(ntracks * sizeof(GXTrack), blockSize);

  // sort SoA tracks by energy in ascending order.
  CopySoATracksToAoS(SoA, AoS, ntracks);
  SortAoSTracksByEnergy(AoS, ntracks);
  CopyAoSTracksToSoA(AoS, SoA, ntracks);

  _mm_free(AoS);
}

void GXTrackHandler::CopyAoSTracks(GXTrack *fromAoS, GXTrack *toAoS,
                                   size_t ntracks)
{
  for (size_t i = 0; i < ntracks; ++i) {
    toAoS[i].x = fromAoS[i].x;
    toAoS[i].y = fromAoS[i].y;
    toAoS[i].z = fromAoS[i].z;
    toAoS[i].t = fromAoS[i].t;
    toAoS[i].px = fromAoS[i].px;
    toAoS[i].py = fromAoS[i].py;
    toAoS[i].pz = fromAoS[i].pz;
    toAoS[i].E = fromAoS[i].E;
    toAoS[i].q = fromAoS[i].q;
    toAoS[i].m = fromAoS[i].m;
  }
}

void GXTrackHandler::CopySoATracks(GXTrack_v &fromSoA, GXTrack_v &toSoA, 
                                   size_t ntracks)
{
  toSoA.size = fromSoA.size;
  for (size_t i = 0; i < ntracks; ++i) {
    (toSoA.x)[i] = (fromSoA.x)[i];
    (toSoA.y)[i] = (fromSoA.y)[i];
    (toSoA.z)[i] = (fromSoA.z)[i];
    (toSoA.t)[i] = (fromSoA.t)[i];
    (toSoA.px)[i] = (fromSoA.px)[i];
    (toSoA.py)[i] = (fromSoA.py)[i];
    (toSoA.pz)[i] = (fromSoA.pz)[i];
    (toSoA.E)[i] = (fromSoA.E)[i];
    (toSoA.q)[i] = (fromSoA.q)[i];
    (toSoA.m)[i] = (fromSoA.m)[i];
  }
}

void GXTrackHandler::CopyAoSTracksToSoA(GXTrack *fromAoS, GXTrack_v &toSoA,
                                        size_t ntracks)
{
  toSoA.size = ntracks;
  for (size_t i = 0; i < ntracks; ++i) {
    (toSoA.x)[i] = fromAoS[i].x;
    (toSoA.y)[i] = fromAoS[i].y;
    (toSoA.z)[i] = fromAoS[i].z;
    (toSoA.t)[i] = fromAoS[i].t;
    (toSoA.px)[i] = fromAoS[i].px;
    (toSoA.py)[i] = fromAoS[i].py;
    (toSoA.pz)[i] = fromAoS[i].pz;
    (toSoA.E)[i] = fromAoS[i].E;
    (toSoA.q)[i] = fromAoS[i].q;
    (toSoA.m)[i] = fromAoS[i].m;
  }
}

void GXTrackHandler::CopySoATracksToAoS(GXTrack_v &fromSoA, GXTrack *toAoS, 
                                        size_t ntracks)
{
  for (size_t i = 0; i < ntracks; ++i) {
    toAoS[i].x = (fromSoA.x)[i];
    toAoS[i].y = (fromSoA.y)[i];
    toAoS[i].z = (fromSoA.z)[i];
    toAoS[i].t = (fromSoA.t)[i];
    toAoS[i].px = (fromSoA.px)[i];
    toAoS[i].py = (fromSoA.py)[i];
    toAoS[i].pz = (fromSoA.pz)[i];
    toAoS[i].E = (fromSoA.E)[i];
    toAoS[i].q = (fromSoA.q)[i];
    toAoS[i].m = (fromSoA.m)[i];
  }
}

} // end namespace gxbert
