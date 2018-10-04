//
// @File: GXCollisionOutput.h
//
// 20180921  Guilherme Lima -- Created, based on M.Kelsey's G4CollisionOutput

#ifndef GXCOLLISION_OUTPUT_HH
#define GXCOLLISION_OUTPUT_HH

//#include "GXFragment.hh"
#include "GXInuclElementaryParticle.hh"
//#include "GXInuclNuclei.hh"
//#include "LorentzRotation.hh"
#include "GXLorentzConvertor.hh"
// //#include "G4ReactionProductVector.hh"

//#include <iosfwd>
#include <algorithm>
#include <vector>

class G4CascadParticle;
class GXLorentzConvertor;

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

template <typename T>
class GXCollisionOutput {

private: 
  int verboseLevel;

  std::vector<GXInuclElementaryParticle<T>> outgoingParticles;
  //std::vector<GXInuclNuclei<T>> outgoingNuclei;
  //std::vector<GXFragment<T>> recoilFragments;
  static const GXFragment emptyFragment;	// To return if list empty

  std::pair<std::pair<int,int>, int> selectPairToTune(double de) const; 
  bool tuneSelectedPair(LorentzVector<T>& mom1, LorentzVector<T>& mom2, int mom_index) const;

  double eex_rest;		// Used by setOnShell() for kinematics
  LorentzVector<T> mom_non_cons;
  bool on_shell;

public:
  GXCollisionOutput() {}
  ~GXCollisionOutput() {}

  GXCollisionOutput& operator=(GXCollisionOutput const& right);

  void setVerboseLevel(int verbose) { verboseLevel = verbose; };

  // ===== Accumulate contents of lists =====

  // Empties lists for new event
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void reset() {
    //outgoingNuclei.clear();
    outgoingParticles.clear();
    //recoilFragments.clear();
    eex_rest = 0.;
    on_shell = false;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void add(const GXCollisionOutput& right);	// Merge complete objects

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void addOutgoingParticle(GXInuclElementaryParticle<T> const& particle)
  {
    outgoingParticles.push_back(particle);
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void addOutgoingParticles(std::vector<GXInuclElementaryParticle<T>> const& particles)
  {
    outgoingParticles.insert(outgoingParticles.end(), particles.begin(), particles.end());
  }

  // void addOutgoingNucleus(const G4InuclNuclei& nuclei) {
  //   outgoingNuclei.push_back(nuclei);
  // };

  // void addOutgoingNuclei(const std::vector<G4InuclNuclei>& nuclea);

  // // These are primarily for G4IntraNucleiCascader internal checks
  // void addOutgoingParticle(const G4CascadeParticle& cparticle);
  // void addOutgoingParticles(const std::vector<G4CascadeParticle>& cparticles);

  //TK for GXBERT
  //void addOutgoingParticles(const G4ReactionProductVector* rproducts);

  // // Special buffer for initial, possible unstable fragments from cascade
  // void addRecoilFragment(const GXFragment* aFragment) {
  //   if (aFragment) addRecoilFragment(*aFragment);
  // }

  // void addRecoilFragment(const GXFragment& aFragment) {
  //   recoilFragments.push_back(aFragment);
  // }

  // ===== Remove contents of lists, by index, reference or value  =====

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void removeOutgoingParticle(int index);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void removeOutgoingParticle(GXInuclElementaryParticle<T> const& particle);

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void removeOutgoingParticle(GXInuclElementaryParticle<T> const* particle)
  {
    if (particle) removeOutgoingParticle(*particle);
  }

  // void removeOutgoingNucleus(G4int index);
  // void removeOutgoingNucleus(const G4InuclNuclei& nuclei);
  // void removeOutgoingNucleus(const G4InuclNuclei* nuclei) {
  //   if (nuclei) removeOutgoingNucleus(*nuclei);
  // }

  // void removeRecoilFragment(G4int index=-1);	// No argument removes all

  // ===== Access contents of lists =====

  int numberOfOutgoingParticles() const { return outgoingParticles.size(); }

  std::vector<GXInuclElementaryParticle<T>> const& getOutgoingParticles() const
  {
    return outgoingParticles;
  };

  std::vector<GXInuclElementaryParticle<T>>& getOutgoingParticles() {
    return outgoingParticles;
  };

  G4int numberOfOutgoingNuclei() const { return 0; } //outgoingNuclei.size(); };
 
  // const std::vector<G4InuclNuclei>& getOutgoingNuclei() const {
  //   return outgoingNuclei;
  // };

  // std::vector<GXInuclNuclei>& getOutgoingNuclei() { return outgoingNuclei; };

  // int numberOfFragments() const { return recoilFragments.size(); }

  // const GXFragment& getRecoilFragment(G4int index=0) const;

  // const std::vector<GXFragment>& getRecoilFragments() const {
  //   return recoilFragments;
  // };

  // std::vector<GXFragment>& getRecoilFragments() { return recoilFragments; };

  // ===== Get event totals for conservation checking, recoil, etc. ======

  LorentzVector<T> getTotalOutputMomentum() const;
  int getTotalCharge() const;			// NOTE:  No fractional charges!
  int getTotalBaryonNumber() const;
  int getTotalStrangeness() const;

  void printCollisionOutput(std::ostream& os = std::cerr) const;

  // ===== Manipulate final-state particles for kinematics =====

  void boostToLabFrame(const GXLorentzConvertor<T>& convertor);

  G4LorentzVector boostToLabFrame(LorentzVector<T>& mom,	// Note pass by value!
				  GXLorentzConvertor<T> const& convertor) const;

  //void rotateEvent(GXLorentzRotation<T> const& rotate);
  void trivialise(GXInuclParticle<T>* bullet, GXInuclParticle<T>* target);
  void setOnShell(GXInuclParticle<T>* bullet, GXInuclParticle<T>* target);
  void setRemainingExcitationEnergy();

  double getRemainingExcitationEnergy() const { return eex_rest; };
  G4bool acceptable() const { return on_shell; };
};

}
}
#endif // GXCOLLISION_OUTPUT_HH
