//
// @File: GXCollisionOutput.h
//
// 20180705  Guilherme Lima -- Created, based on M.Kelsey's G4CollisionOutput

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

class GXCollisionOutput {

private: 
  int verboseLevel;

  std::vector<GXInuclElementaryParticle<double>> outgoingParticles;
  //std::vector<GXInuclNuclei<double>> outgoingNuclei;
  //std::vector<GXFragment<double>> recoilFragments;
  static const GXFragment emptyFragment;	// To return if list empty

  std::pair<std::pair<int,int>, int> selectPairToTune(double de) const; 
  bool tuneSelectedPair(LorentzVector<double>& mom1, LorentzVector<double>& mom2,
			int mom_index) const;

  double eex_rest;		// Used by setOnShell() for kinematics
  LorentzVector<double> mom_non_cons;
  bool on_shell;

public:
  GXCollisionOutput() {}
  ~GXCollisionOutput() {}

  GXCollisionOutput& operator=(const GXCollisionOutput& right);

  void setVerboseLevel(G4int verbose) { verboseLevel = verbose; };

  // ===== Accumulate contents of lists =====

  void reset();		// Empties lists for new event

  void add(const GXCollisionOutput& right);	// Merge complete objects

  void addOutgoingParticle(const GXInuclElementaryParticle<double>& particle) {
    outgoingParticles.push_back(particle);
  }

  void addOutgoingParticles(const std::vector<GXInuclElementaryParticle<double>>& particles);

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

  void removeOutgoingParticle(int index);
  void removeOutgoingParticle(const GXInuclElementaryParticle<double>& particle);
  void removeOutgoingParticle(const GXInuclElementaryParticle<double>* particle) {
    if (particle) removeOutgoingParticle(*particle);
  }

  // void removeOutgoingNucleus(G4int index);
  // void removeOutgoingNucleus(const G4InuclNuclei& nuclei);
  // void removeOutgoingNucleus(const G4InuclNuclei* nuclei) {
  //   if (nuclei) removeOutgoingNucleus(*nuclei);
  // }

  // void removeRecoilFragment(G4int index=-1);	// No argument removes all

  // ===== Access contents of lists =====

  G4int numberOfOutgoingParticles() const { return outgoingParticles.size(); }
    
  const std::vector<GXInuclElementaryParticle<double>>& getOutgoingParticles() const {
    return outgoingParticles;
  };

  std::vector<GXInuclElementaryParticle<double>>& getOutgoingParticles() {
    return outgoingParticles;
  };

  //G4int numberOfOutgoingNuclei() const { return outgoingNuclei.size(); };
 
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

  G4LorentzVector getTotalOutputMomentum() const;
  G4int getTotalCharge() const;			// NOTE:  No fractional charges!
  G4int getTotalBaryonNumber() const;
  G4int getTotalStrangeness() const;

  void printCollisionOutput(std::ostream& os=G4cout) const;

  // ===== Manipulate final-state particles for kinematics =====

  void boostToLabFrame(const GXLorentzConvertor<double>& convertor);

  G4LorentzVector boostToLabFrame(LorentzVector<double>& mom,	// Note pass by value!
				  GXLorentzConvertor<double> const& convertor) const;

  //void rotateEvent(GXLorentzRotation<T> const& rotate);
  void trivialise(GXInuclParticle<double>* bullet, GXInuclParticle<double>* target);
  void setOnShell(GXInuclParticle<double>* bullet, GXInuclParticle<double>* target);
  void setRemainingExcitationEnergy();

  double getRemainingExcitationEnergy() const { return eex_rest; };
  G4bool acceptable() const { return on_shell; };
};        

}
}
#endif // GXCOLLISION_OUTPUT_HH
