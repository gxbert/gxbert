//
// @File: GXVTwoBodyAngDst.hh
//
// Description: pure virtual base class for two-body final state angular
//              distributions in Bertini-style cascade
//
// 20130221  M. Kelsey -- Bug fix: theName should not be a reference,
//		but a local copy.
// 20130308  M. Kelsey -- Add access to name string for diagnostic utilities
// 2018-080-15  Guilherme Lima  - Created, based on M.Kelsey's G4VTwoBodyAngDst.hh
//
#ifndef GXVTwoBodyAngDst_h
#define GXVTwoBodyAngDst_h 1

#include <string>

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

template <typename T>
class GXVTwoBodyAngDst {
public:
  GXVTwoBodyAngDst(const std::string& name, int verbose=0);

  virtual ~GXVTwoBodyAngDst() { }

  virtual T GetCosTheta(const T& ekin, const T& pcm) const = 0;

  virtual void setVerboseLevel(int verbose = 0) { verboseLevel = verbose; }

  virtual const std::string& GetName() const { return theName; }

protected:
  std::string theName;
  int verboseLevel;
};        

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXVTwoBodyAngDst_h
