#ifndef StMixedSoftPion_hh
#define StMixedSoftPion_hh
/* **************************************************
 *
 * Track class used for mixed event buffer, stripped down
 * to minimum information neede to reconstruct the helix
 * and basic track information. Currently include:
 * 1) charge
 * 2) isTpcPi & isTofPi
 * 3) isTpcKaon & is TofKaon
 *
 * **************************************************
 *
 *  Initial Authors:
 *         ** Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */
#include <math.h>
#include "TVector2.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"

class StPicoTrack;

class StMixedSoftPion : public TObject
{
public:
   StMixedSoftPion();
   StMixedSoftPion(StMixedSoftPion const *);
   StMixedSoftPion(StThreeVectorF const & pVtx, float B, StPicoTrack const* const picoTrack);
   short const charge() const ;
   StThreeVectorF const& gMom() const;
   StThreeVectorF const& origin() const;
   StLorentzVectorF const& lorentzVector() const;
   float const pt() const { return mPt; }
   ~StMixedSoftPion()
   {
      ;
   };
private:
   StThreeVectorF mOrigin;
   StThreeVectorF mMom;
   short mCharge;
   StLorentzVectorF mLorentzVector;
   float mPt;
   //Removed origin, allt racks shoud me set to 0,0,0
   ClassDef(StMixedSoftPion,1);
};

inline StThreeVectorF const & StMixedSoftPion::gMom() const
{
   return (mMom) ;
}
inline StThreeVectorF const & StMixedSoftPion::origin() const
{
   return (mOrigin) ;
}
inline short const StMixedSoftPion::charge() const
{
   return mCharge;
}
inline StLorentzVectorF const& StMixedSoftPion::lorentzVector() const
{
   return mLorentzVector;
}
#endif
