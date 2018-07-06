#ifndef GXINUCL_SPECIAL_FUNCTIONS_HH
#define GXINUCL_SPECIAL_FUNCTIONS_HH 1

/* 
SIMD/SIMT version of G4InuclSpecialFunctions.hh/cc
*/

#include "VecHepDefs.hh"
#include "GXThreeVector.hh"
#include "LorentzVector.hh"
#include "GXRandom.hh"
#include "CLHEP/Units/SystemOfUnits.h"

namespace gxbert {
inline namespace GXBERT_IMPL_NAMESPACE {

namespace GXInuclSpecialFunctions {

  // csNN
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v csNN(typename Backend::Double_v e) 
  {

    using Double_v = typename VectorBackend::Double_v;

    Mask_v<Double_v>  low = (e < 40.0);
    Double_v snn = Blend(low, -1174.8 / (e * e) + 3088.5 / e + 5.3107,
                              93074.0 / (e * e) - 11.148 / e + 22.429);
    return snn;
  }

  // csPN
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v csPN(typename Backend::Double_v e) 
  {
    using Double_v = typename VectorBackend::Double_v;

    Mask_v<Double_v>  low = (e < 40.0);
    Double_v spn = Blend(low,  -5057.4 / (e * e) + 9069.2 / e + 6.9466,
                              239380.0 / (e * e) + 1802.0 / e + 27.147);
    return spn;
  }

  // FermiEnergy
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v FermiEnergy(Index_v<typename Backend::Double_v> A,
                                         Index_v<typename Backend::Double_v> Z,
                                         Index_v<typename Backend::Double_v> ntype)
  {
    using Double_v = typename VectorBackend::Double_v;

    Double_v Z23 = 2./3.; 
    Double_v C =  55.4 / math::Pow(math::Log(Double_v(math::Abs(A))), Z23);

    Mask_v<Index_v<Double_v>>  nzero = (ntype == 0);
    Double_v arg  = Blend(nzero, math::Pow(math::Log(Double_v(math::Abs(A-Z))),Z23),
			  math::Pow(math::Log(Double_v(math::Abs(Z))),Z23));
    return C * arg;
  }

  // GXcbrt
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v GXcbrt(typename Backend::Double_v x) 
  {
    using Double_v = typename Backend::Double_v;
    Mask_v<Double_v>  nonzero = (x != 0.);
    Mask_v<Double_v>  negative = (x < 0.);
    Double_v fact  = Blend(nonzero,Blend(negative, -1., 1.),0.);
    return fact*math::Exp(math::Log(math::Abs(x))/3.);
  }

  // GXcbrt - int
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v GXcbrt(Index_v<typename Backend::Double_v> n) 
  {
    using Double_v = typename Backend::Double_v;
    Mask_v<Index_v<Double_v>>  nonzero = (n != 0.);
    Mask_v<Index_v<Double_v>>  negative = (n < 0.);
    Double_v fact  = Blend(nonzero,Blend(negative, -1., 1.),0.);
    return fact*math::Pow(math::Log(math::Abs(n)),1/3.);
  }

  // getAL
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v getAL(Index_v<typename Backend::Double_v> A)
  {
    return 0.76 + 2.2 / GXcbrt<Backend>(A);
  }

  // inuclRndm
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v inuclRndm() 
  { 
    return gxbert::cxx::GXRandom::RNG().Uniform<Backend>();
  } 

  // randomPHI
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  typename Backend::Double_v randomPHI() 
  { 
    return CLHEP::twopi * inuclRndm<Backend>();
  } 

  // randomCOS_SIN
  template <class Backend>
  std::pair<typename Backend::Double_v, typename Backend::Double_v> randomCOS_SIN() 
  {
    using Double_v = typename VectorBackend::Double_v;
    Double_v CT = 1.0 - 2.0 * inuclRndm<Backend>();
    return std::pair<Double_v, Double_v>(CT, math::Sqrt(1.0 - CT*CT));
  }

  // randomGauss
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  typename Backend::Double_v randomGauss(typename Backend::Double_v sigma) 
  {
     /*
     const G4double eps = 1.0e-6;
     G4double r1 = inuclRndm();
     r1 = r1 > eps ? r1 : eps;
     G4double r2 = inuclRndm();
     r2 = r2 > eps ? r2 : eps;
     r2 = r2 < 1.0 - eps ? r2 : 1.0 - eps; 
 
     return sigma * std::sin(twopi * r1) * std::sqrt(-2.0 * GXLog(r2)); 
     */

    //use Gauss variate directly for the generator
    return gxbert::cxx::GXRandom::RNG().Gauss<Backend>(0.0,sigma);
  } 

  // generateWithFixedTheta
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  LorentzVector<typename Backend::Double_v> generateWithFixedTheta(typename Backend::Double_v ct, 
                                                                   typename Backend::Double_v p, 
                                                                   typename Backend::Double_v mass) 
  {
    using Double_v = typename Backend::Double_v;

    Double_v phi = randomPHI<Backend>();
    Double_v pt = p * math::Sqrt(math::Abs(1.0 - ct * ct));

    GXThreeVector<Double_v> pvec;
    LorentzVector<Double_v> momr;

    pvec.Set(pt*math::Cos(phi), pt*math::Sin(phi), p*ct);
    momr.SetVectM(pvec, mass);

    return momr;
  }

  // generateWithRandomAngles
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  LorentzVector<typename Backend::Double_v> generateWithRandomAngles(typename Backend::Double_v p, 
                                                                     typename Backend::Double_v mass) 
  {
    using Double_v = typename Backend::Double_v;

    std::pair<Double_v, Double_v> COS_SIN = randomCOS_SIN<Backend>();

    Double_v phi = randomPHI<Backend>();
    Double_v pt = p * COS_SIN.second;

    GXThreeVector<Double_v> pvec;
    LorentzVector<Double_v> momr;

    pvec.Set(pt*math::Cos(phi), pt*math::Sin(phi), p*COS_SIN.first);
    momr.SetVectM(pvec, mass);

    return momr;
  }

} // end namespace G4InuclSpecialFunctions
} // end namespace impl
} // end namespace vecCore

#endif
