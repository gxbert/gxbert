#ifndef GXTRACK_H
#define GXTRACK_H 1

#include "stddef.h"
#include "GXThreeVector.hh"
#include "LorentzVector.hh"
#include "G4LorentzVector.hh"

/**
 * GXTrack: A utility struct of input/output data for unit/benchmark tests
 */

namespace gxbert {

struct GXTrack {
  double x;     // position x
  double y;     // position y
  double z;     // position z
  double t;     // time
  double px;    // momentum x
  double py;    // momentum y
  double pz;    // momentum z
  double E;     // total energy
  double q;     // charge
  double m;     // mass
};

//SOA 

struct GXTrack_v {
  int    size;  // number of tracks
  double *x;
  double *y;
  double *z;
  double *t;
  double *px;
  double *py;
  double *pz;
  double *E; 
  double *q;
  double *m;

  template <typename Real_v>
  void getThreeMomentum(size_t i, GXThreeVector<Real_v> &out) const { out.Set(*(Real_v*)(px + i), *(Real_v*)(py + i), *(Real_v*)(pz + i)); }

  template <typename Real_v>
  GXThreeVector<Real_v> getThreeMomentum(size_t i) const { return GXThreeVector<Real_v>(*(Real_v*)(px+i), *(Real_v*)(py+i), *(Real_v*)(pz+i)); }

  template <typename Real_v>
  void getFourMomentum(size_t i, LorentzVector<Real_v> &out) const { out.Set(*(Real_v*)(px + i), *(Real_v*)(py + i), *(Real_v*)(pz + i), *(Real_v*)(E + i) + *(Real_v*)(m + i)); }

  template <typename Real_v>
  LorentzVector<Real_v> getFourMomentum(size_t i) const { return LorentzVector<Real_v>(*(Real_v*)(px+i), *(Real_v*)(py+i), *(Real_v*)(pz+i), *(Real_v*)(E+i) + *(Real_v*)(m+i)); }

  // for CLHEP compatibility
  void getFourMomentum(size_t i, G4LorentzVector &out) const
  {
    out.setX(*(px + i));
    out.setY(*(py + i));
    out.setZ(*(pz + i));
    out.setT(*(E + i) + *(m + i));
  }

};

} // end namespace gxbert

#endif
