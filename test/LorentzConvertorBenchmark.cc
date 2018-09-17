
#include "VecCore/VecCore"
#include "timer.h"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzConvertor.hh"

#include "GXTrack.hh"
#include "GXTrackHandler.hh"
#include "GXLorentzConvertor.hh"
#include "GXProton.hh"
#include "GXNeutron.hh"
#include "LorentzVector.hh"

using namespace gxbert;

//.. default globals
constexpr double pmass = 938.272013;  // proton mass in MeV
int debugLevel   = 0;
int nReps        = 1;
int nEvents      = 1024*64;
double kinEnergy = 1.3E+3; // in MeV

void RunG4LorentzConvertorTimer(GXTrack_v const& soaBullets, GXTrack_v const& soaTargets)
{
  //==== Reference: G4LorentzConvertor functionalities
  double g4sumEscm = 0.;
  double g4sumEkin = 0.;
  double g4sumP2 = 0.;
  G4LorentzConvertor g4conv;
  G4LorentzVector lorvec;
  g4conv.setVerbose(debugLevel);

  Timer<nanoseconds> timer;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    g4sumEscm = g4sumEkin = g4sumP2 = 0.;
    for(int i = 0; i < nEvents; ++i) {
      soaTargets.getFourMomentum(i, lorvec);
      g4conv.setTarget(lorvec);

      soaBullets.getFourMomentum(i, lorvec);
      g4conv.setBullet(lorvec);

      g4conv.toTheCenterOfMass();
      g4sumEscm += g4conv.getTotalSCMEnergy();
      g4sumEkin += g4conv.getKinEnergyInTheTRS();

      G4LorentzVector plab = g4conv.backToTheLab(lorvec);
      g4sumP2 += plab.mag2();
    }
  }
  double g4Elapsed = timer.Elapsed();

  std::cerr<<"GXBert results: "
	   <<"  sumEscm = "<< g4sumEscm
	   <<"\t  sumEkin = "<< g4sumEkin
	   <<"\t  sumP2 = "<< g4sumP2
	   <<"\t CPUtime = "<< g4Elapsed / nEvents / nReps
	   <<"\n";
}


template <typename Real_v>
void RunGXLorentzConvertorTimer(const char* testname, GXTrack_v const& soaBullets, GXTrack_v const& soaTargets)
{
  int vsize = vecCore::VectorSize<Real_v>();
  Real_v sumEscm = 0.;
  Real_v sumEkin = 0.;
  Real_v sumP2 = 0.;
  GXLorentzConvertor<Real_v> conv;
  conv.setVerbose(debugLevel);

  LorentzVector<Real_v> fourvec;
  Timer<nanoseconds> timer;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    sumEscm = sumEkin = sumP2 = 0.;
    for(int i = 0; i < nEvents; i += vsize) {
      soaTargets.getFourMomentum(i, fourvec);
      conv.setTarget(fourvec);

      soaBullets.getFourMomentum(i, fourvec);
      conv.setBullet(fourvec);

      conv.toTheCenterOfMass();
      sumEscm += conv.getTotalSCMEnergy();
      sumEkin += conv.getKinEnergyInTheTRS();

      LorentzVector<Real_v> plab = conv.backToTheLab(fourvec);
      sumP2 += plab.Mag2();
    }
  }
  double elapsed = timer.Elapsed();

  std::cerr<< testname <<" results: "
	   <<   "sumEscm = "<< vecCore::ReduceAdd(sumEscm)
	   <<"\t sumEkin = "<< vecCore::ReduceAdd(sumEkin)
	   <<"\t sumP2 = "  << vecCore::ReduceAdd(sumP2)
	   <<"\t CPUtime = "<< elapsed / nEvents / nReps
	   <<"\n";
}


int main(int argc, char* argv[]) {

  if(argc >= 2) kinEnergy  = atof(argv[1]);
  if(argc >= 3) nEvents = atoi(argv[2]);
  if(argc >= 4) nReps = atoi(argv[3]);

  const double etot = kinEnergy + pmass;
  double pmom = sqrt(etot * etot - pmass * pmass);

  //.. timer
  Timer<nanoseconds> timer;

  //.. prepare bullets
  GXTrackHandler *bulletHandler = new GXTrackHandler(nEvents);
  bulletHandler->GenerateRandomTracks(nEvents, pmom, pmom);
  GXTrack_v &soaBullets = bulletHandler->GetSoATracks();
  //printf("Number of tracks = %d\n", soaBullets.size);

  //.. prepare targets
  GXTrackHandler *targetHandler = new GXTrackHandler(nEvents);
  targetHandler->GenerateRandomTracks(nEvents, pmom, pmom);
  GXTrack_v &soaTargets = targetHandler->GetSoATracks();

  //==== Reference: G4LorentzConvertor functionalities
  double g4sumEscm = 0.;
  double g4sumEkin = 0.;
  double g4sumP2 = 0.;
  G4LorentzConvertor g4conv;
  G4LorentzVector lorvec;
  g4conv.setVerbose(debugLevel);

  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    g4sumEscm = g4sumEkin = g4sumP2 = 0.;
    for(int i = 0; i < nEvents; ++i) {
      soaTargets.getFourMomentum(i, lorvec);
      g4conv.setTarget(lorvec);

      soaBullets.getFourMomentum(i, lorvec);
      g4conv.setBullet(lorvec);

      g4conv.toTheCenterOfMass();
      g4sumEscm += g4conv.getTotalSCMEnergy();
      g4sumEkin += g4conv.getKinEnergyInTheTRS();

      G4LorentzVector plab = g4conv.backToTheLab(lorvec);
      g4sumP2 += plab.mag2();
    }
  }
  double g4Elapsed = timer.Elapsed();

  std::cerr<<"GXBert results: "
	   <<"  sumEscm = "<< g4sumEscm
	   <<"\t  sumEkin = "<< g4sumEkin
	   <<"\t  sumP2 = "<< g4sumP2
	   <<"\t CPUtime = "<< g4Elapsed / nEvents / nReps
	   <<"\n";

  //=== Templated converter in scalar mode
  double gsSumEscm = 0.;
  double gsSumEkin = 0.;
  double gsSumP2 = 0.;
  GXLorentzConvertor<double> gsconv;
  gsconv.setVerbose(debugLevel);

  LorentzVector<double> fourvec;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    gsSumEscm = gsSumEkin = gsSumP2 = 0.;
    for(int i = 0; i < nEvents; i++) {
      soaTargets.getFourMomentum(i, fourvec);
      gsconv.setTarget(fourvec);

      soaBullets.getFourMomentum(i, fourvec);
      gsconv.setBullet(fourvec);

      gsconv.toTheCenterOfMass();
      gsSumEscm += gsconv.getTotalSCMEnergy();
      gsSumEkin += gsconv.getKinEnergyInTheTRS();

      LorentzVector<double> plab = gsconv.backToTheLab(fourvec);
      gsSumP2 += plab.Mag2();
    }
  }
  double gsElapsed = timer.Elapsed();

  std::cerr<<"Scalar results: "
	   <<"  sumEscm = "<< gsSumEscm
	   <<"\t  sumEkin = "<< gsSumEkin
	   <<"\t  sumP2 = "<< gsSumP2
	   <<"\t CPUtime = "<< gsElapsed / nEvents / nReps
	   <<"\n";


  //=== Vectorized converter
  //.. vector size 
  using Real_v = typename vecCore::backend::VcVector::Double_v;
  int vsize = vecCore::VectorSize<Real_v>();
  std::cerr<<"Vector size: "<< vsize <<"\n";

  Real_v gxSumEscm = 0.;
  Real_v gxSumEkin = 0.;
  Real_v gxSumP2 = 0.;
  GXLorentzConvertor<Real_v> gxconv;
  gxconv.setVerbose(debugLevel);

  LorentzVector<Real_v> simd4vec;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    gxSumEscm = gxSumEkin = gxSumP2 = 0.;
    for(int i = 0; i < nEvents; i += vsize) {
      soaTargets.getFourMomentum(i, simd4vec);
      gxconv.setTarget(simd4vec);

      soaBullets.getFourMomentum(i, simd4vec);
      gxconv.setBullet(simd4vec);

      gxconv.toTheCenterOfMass();
      gxSumEscm += gxconv.getTotalSCMEnergy();
      gxSumEkin += gxconv.getKinEnergyInTheTRS();

      LorentzVector<Real_v> plab = gxconv.backToTheLab(simd4vec);
      gxSumP2 += plab.Mag2();
    }
  }
  double gxElapsed = timer.Elapsed();

  std::cerr<<"Vector results: "
	   <<"  sumEscm = "<< vecCore::ReduceAdd(gxSumEscm)
	   <<"\t  sumEkin = "<< vecCore::ReduceAdd(gxSumEkin)
	   <<"\t  sumP2 = "<< vecCore::ReduceAdd(gxSumP2)
	   <<"\t CPUtime = "<< gxElapsed / nEvents / nReps
	   <<"\n";


  //=== Alternate way to load arrays for vectorized converter
  Double_v px, py, pz, m, kinE;
  Real_v glSumEscm = 0.;
  Real_v glSumEkin = 0.;
  Real_v glSumP2 = 0.;
  timer.Start();
  for(int irep = 0; irep < nReps; ++irep) {
    glSumEscm = glSumEkin = glSumP2 = 0.;
    for (int i = 0; i < nEvents; i += vsize) {
      Load(px, &soaTargets.px[i]);
      Load(py, &soaTargets.py[i]);
      Load(pz, &soaTargets.pz[i]);
      Load( m, &soaTargets.m[i]);
      Load(kinE, &soaTargets.E[i]);
      simd4vec.Set(px, py, pz, kinE+m);
      gxconv.setTarget(simd4vec);//, simd4vec.Mag());

      Load(px, &soaBullets.px[i]);
      Load(py, &soaBullets.py[i]);
      Load(pz, &soaBullets.pz[i]);
      Load( m, &soaBullets.m[i]);
      Load(kinE, &soaBullets.E[i]);
      simd4vec.Set(px, py, pz, kinE+m);
      gxconv.setBullet(simd4vec);//, simd4vec.Mag());

      gxconv.toTheCenterOfMass();
      glSumEscm += gxconv.getTotalSCMEnergy();
      glSumEkin += gxconv.getKinEnergyInTheTRS();

      LorentzVector<Real_v> plab = gxconv.backToTheLab(simd4vec);
      glSumP2 += plab.Mag2();
    }
  }
  double glElapsed = timer.Elapsed();

  std::cerr<<"VectorL result: "
	   <<"  sumEscm = "<< vecCore::ReduceAdd(glSumEscm)
	   <<"\t  sumEkin = "<< vecCore::ReduceAdd(glSumEkin)
	   <<"\t  sumP2 = "<< vecCore::ReduceAdd(glSumP2)
	   <<"\t CPUtime = "<< glElapsed / nEvents / nReps
	   <<"\n\n\n";

  // Using function calls for benchmarks
  RunG4LorentzConvertorTimer(soaBullets, soaTargets);

  // benchmarks templated on types
  RunGXLorentzConvertorTimer<double>("double", soaBullets, soaTargets);
  RunGXLorentzConvertorTimer<vecCore::backend::VcVector::Double_v>("Double_v", soaBullets, soaTargets);

  //RunGXLorentzConvertorTimer<float>("float", soaBullets, soaTargets);
  //RunGXLorentzConvertorTimer<vecCore::backend::VcVector::Float_v>("Double_v", soaBullets, soaTargets);

  //.. cleanup
  delete bulletHandler;
  delete targetHandler;
}
