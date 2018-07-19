//
//
#ifndef GXCascadeFinalStateGenerator_hh
#define GXCascadeFinalStateGenerator_hh 1

//#include "globals.hh"
#include "GXCascadeGenerator.hh"
#include "GXCascadeFinalStateAlgorithm.hh"
#include <vector>

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

  template <typename T>
  class GXInuclElementaryParticle;

  template <typename T>
  class GXCascadeFinalStateAlgorithm;

  //=== class definition ===
  template <typename T>
  class GXCascadeFinalStateGenerator : public GXCascadeGenerator<T> {
  public:
    GXCascadeFinalStateGenerator()
      : GXCascadeGenerator<T>(new GXCascadeFinalStateAlgorithm<T>)
    { }

    virtual ~GXCascadeFinalStateGenerator() {}

    //.. Configure algorithm (distributions) based on interaction
    void Configure(GXInuclElementaryParticle<T> const* bullet,
		   GXInuclElementaryParticle<T> const* target,
		   const std::vector<Index_v<T>>& particle_kinds);
  };

  // === member declarations
  template <typename T>
  void GXCascadeFinalStateGenerator<T>::Configure(GXInuclElementaryParticle<T> const* bullet,
						  GXInuclElementaryParticle<T> const* target,
						  const std::vector<Index_v<T>>& particle_kinds)
  {
    return;
  }
} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif	// GXCascadeFinalStateGenerator_hh
