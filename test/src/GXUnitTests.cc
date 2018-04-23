template <typename T>                                       
VECCORE_ATT_HOST_DEVICE 
VECCORE_FORCE_INLINE
LorentzVector<T>& LorentzVector<T>::Boost(T bx, T by, T bz)
{
  T b2 = bx*bx + by*by + bz*bz;
  T ggamma = 1.0 / math::Sqrt(1.0 - b2);
  T bp = bx*x() + by*y() + bz*z();

  Mask_v<T> positive = (b2 > 0.);
  T gamma2 =  Blend(positive,(ggamma - 1.0)/b2, static_cast<T>(0.0));

  SetX(x() + gamma2*bp*bx + ggamma*bx*t());
  SetY(y() + gamma2*bp*by + ggamma*by*t());
  SetZ(z() + gamma2*bp*bz + ggamma*bz*t());
  SetT(ggamma*(t() + bp));
  return *this;
} 
