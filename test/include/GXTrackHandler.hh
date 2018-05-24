#ifndef GXTrackHandler_H
#define GXTrackHandler_H 1

#include "GXTrack.hh"
#include "stddef.h"

#include "VecRng/MRG32k3a.h"
#include "VecHepDefs.hh"

namespace gxbert {

class GXTrackHandler {
public:
  GXTrackHandler();
  GXTrackHandler(size_t nTracks);
  ~GXTrackHandler();

  void SetNumberOfTracks(size_t nTracks); 
  size_t GetNumberOfTracks() { return fNumberOfTracks; };

  void Allocate(size_t nTracks);
  void Deallocate();
  void Reallocate(size_t nTracks);
  void SetMass(double mass) { fMass = mass; }

  double Random();
  void SetRandomStream(long streamId);

  GXTrack *GetAoSTracks() { return fTrack_aos; };
  GXTrack_v &GetSoATracks() { return fTrack_soa; };

  void GenerateRandomTracks(size_t nTracks, double minP = 1., double maxP = 1E+4);

  // generate nTracks at a given position + direction given by posdir[6], with 3-mom/GeV in range (minP, maxP)
  void GenerateTracksAlongSameDirection(size_t nTracks, double* posdir, double minP = 1., double maxP = 1E+4);

  // utility functions
  void SortAoSTracksByEnergy(GXTrack *AoS, size_t nTracks);
  void SortSoATracksByEnergy(GXTrack_v &SoA, size_t nTracks);

  void CopyAoSTracks(GXTrack *fromAoS, GXTrack *toAoS, size_t Tracks);
  void CopySoATracks(GXTrack_v &fromSoA, GXTrack_v &toSoA, size_t nTracks);
  void CopyAoSTracksToSoA(GXTrack *fromAoS, GXTrack_v &toSoA, size_t nTracks);
  void CopySoATracksToAoS(GXTrack_v &fromSoA, GXTrack *toAoS, size_t nTracks);

private:
  size_t fNumberOfTracks;
  GXTrack *fTrack_aos;
  GXTrack_v fTrack_soa;
  char *fBuffer;
  double fMass; // mass to be used for tracks generated (all tracks of the same type)

  vecRng::MRG32k3a<ScalarBackend> *fRNG;
  
};

} // end namespace gxbert

#endif
