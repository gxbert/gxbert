#ifndef GXCURAND_INIT_H
#define GXCURAND_INIT_H 1

#include "VecRng/RngDefs.h"

namespace gxbert {

bool GXCurand_Init(Random_t *randomStates, unsigned long seed, int blocksPerGrid, int threadsPerBlock);

} // end namespace gxbert

#endif
