#ifndef GXTRACK_H
#define GXTRACK_H 1

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
};

} // end namespace gxbert

#endif
