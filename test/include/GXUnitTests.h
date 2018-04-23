#ifndef GXUNITTESTS_H
#define GXUNITTESTS_H 1

#include "VecHepDefs.h"

namespace gxbert {
inline namespace GXBERT_IMPL_NAMESPACE {

template <typename BackendT>
class GXUnitTests 
{

public:

  VECCORE_ATT_HOST_DEVICE
  GXUnitTests() {}

  VECCORE_ATT_HOST_DEVICE
  ~GXUnitTests() {}

  // Exponential deviates: exp(-x/tau)
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v LorentzBoost(typename BackendT::Double_v x,
                                           typename BackendT::Double_v z,
                                           typename BackendT::Double_v y);
};  

} // end namespace impl
} // end namespace gxbert

#endif
