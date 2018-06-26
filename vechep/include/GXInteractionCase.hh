//
// File:    GXInteractionCase.hh
//
// Purpose: Store input objects to Bertini algorithms
//
// 20180606  Guilherme Lima    Created, based on M. Kelsey's G4InteractionCase
//

#ifndef GXInteractionCase_HH
#define GXInteractionCase_HH

#include "globals.hh"
#include "GXInuclParticle.hh"
#include "GXInuclElementaryParticle.hh"
#include <sstream>

namespace gxbert {

inline namespace GXBERT_IMPL_NAMESPACE {

template <typename T>
class GXInteractionCase {

  using Bool_v = vecCore::Mask_v<T>;

private:
  GXInuclParticle<T> const* bullet;
  GXInuclParticle<T> const* target;
  Index_v<T> fCase;

public:
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInteractionCase() : bullet(0), target(0), fCase(0) {}

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInteractionCase(GXInuclParticle<T> const* part1, GXInuclParticle<T> const* part2) {
    set(part1, part2);
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void set(GXInuclParticle<T> const* part1, GXInuclParticle<T> const* part2)
  {
    clear();  // Reset everything in case of failure

    // See which one of the two (or both) is a nucleus
    //GXInuclNuclei<T> const* nucl1 = dynamic_cast<GXInuclNuclei const*>(part1);
    //GXInuclNuclei<T> const* nucl2 = dynamic_cast<GXInuclNuclei const*>(part2);

    GXInuclElementaryParticle<T> const* had1 = dynamic_cast<GXInuclElementaryParticle<T> const*>(part1);
    GXInuclElementaryParticle<T> const* had2 = dynamic_cast<GXInuclElementaryParticle<T> const*>(part2);

    // if (nucl1 && nucl2) { 	// Nuclear collision, lighter is projectile
    //   fCase = -2;
    //   if (nucl2->getA() >= nucl1->getA()) {
    // 	bullet = part1;
    // 	target = part2;
    //   } else {
    // 	bullet = part2;
    // 	target = part1;
    //   } 
    // } else if (nucl1 || nucl2) {	// Hadron on nucleus, hadron projectile
    //   fCase = -1;
    //   if (nuclq1 && had2) {
    // 	bullet q= part2;
    // 	target = part1;
    //   } else {
    // 	bullet = part1;
    // 	target = part2;
    //   }
    //} elseif (had1 && had2) {	// Hadron-hadron interaction, order irrelevant
    if (had1 && had2) {	// Hadron-hadron interaction, order irrelevant
      fCase = had1->type() * had2->type();
      bullet = part1;
      target = part2;
    }
    //std::cerr<<"*** GXInterCase: bullet: "<< *bullet <<"\n";
    //std::cerr<<"*** GXInterCase: target: "<< *target <<"\n";
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void clear() {
    bullet = target = 0;
    fCase = 0;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle<T> const* getBullet() const { return bullet; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  GXInuclParticle<T> const* getTarget() const { return target; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v match(int ityp) const
  {
    Index_v<T> ref( (Index<T>)ityp );
    return fCase == ref;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v valid() const      { return ! match(0); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v twoNuclei() const  { return match(-2); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Bool_v hadNucleus() const { return match(-1); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  Index_v<T> hadrons() const    { return fCase; }	// "rtype" or "is" code

  // // For compatibility with GXIntraNucleiCascader code
  // VECCORE_ATT_HOST_DEVICE
  // VECCORE_FORCE_INLINE
  // Index_v<T>  code() const { return ((fCase<0) ? -fCase : 0); }
};



//=== Debugging helpers
template <typename T>
std::ostream &operator<<(std::ostream &os, typename gxbert::GXInteractionCase<T> const& intcase)
{
  auto vsize = vecCore::VectorSize<T>();
  std::stringstream tnames, bnames;
  std::stringstream tmasses, bmasses;

  GXInuclElementaryParticle<T> const* bhad = dynamic_cast<GXInuclElementaryParticle<T> const*>(intcase.getBullet());
  GXInuclElementaryParticle<T> const* thad = dynamic_cast<GXInuclElementaryParticle<T> const*>(intcase.getTarget()); 
  Index_v<T> btyp = bhad ? bhad->type() : 0;
  Index_v<T> ttyp = thad ? thad->type() : 0;
  for (int i = 0; i < vsize; ++i) {
    if (i > 0) {
      bnames <<"; ";
      bmasses <<"; ";
      tnames <<"; ";
      tmasses <<"; ";
    }
    GXParticleDefinition const* pd = getDefinition(Get(btyp, i));
    bnames << pd->GetParticleName();
    bmasses << pd->GetPDGMass();

    pd = getDefinition(Get(ttyp, i));
    tnames << pd->GetParticleName();
    tmasses << pd->GetPDGMass();
  }
  os << " InterCase bullets=["<< bnames.str() <<"]"
     << " masses=[" << bmasses.str() <<"]"
     << " types=" << btyp
     << " ekin=" << bhad->getKineticEnergy();

  os << " InterCase targets=["<< tnames.str() <<"]"
     << " masses=[" << tmasses.str() <<"]"
     << " types=" << ttyp
     << " ekin=" << thad->getKineticEnergy();

  return os;
}

template <>
std::ostream &operator<<(std::ostream &os, GXInteractionCase<double> const& intcase)
{
  GXInuclElementaryParticle<double> const* trk = dynamic_cast<GXInuclElementaryParticle<double> const*>(intcase.getBullet());
  GXParticleDefinition const* pd = getDefinition(trk->type());
  os << "\n\t=== InterCase bullet: " << pd->GetParticleName()
     << " mass=" << pd->GetPDGMass()
     << " type=" << trk->type()
     << " ekin=" << trk->getKineticEnergy();

  trk = dynamic_cast<GXInuclElementaryParticle <double>const*>(intcase.getTarget());
  pd = getDefinition(trk->type());
  os << "\n\t=== InterCase target: " << pd->GetParticleName()
     << " mass=" << pd->GetPDGMass()
     << " type=" << trk->type()
     << " ekin=" << trk->getKineticEnergy();
  os <<"\n";
  return os;
}

} // end GXBERT_IMPL_NAMESPACE
} // end namespace gxbert
#endif // GXINTERACTIONCASE_HH 
