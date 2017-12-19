#ifndef StMixedD0_hh
#define StMixedD0_hh
/* **************************************************
 *
 * Track class used for mixed event buffer, stripped down
 * to minimum information neede to reconstruct the helix
 * and basic track information. Currently include:
 * 1) charge
 * 2) 
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
#include "StPicoCharmContainers/StKaonPion.h"
class StPicoTrack;
class StMixedD0 : public TObject
{
public:
   StMixedD0();
   StMixedD0(StMixedD0 const *);
   StMixedD0(StThreeVectorF const & pVtx, float B, StKaonPion const * const aD0,StPicoTrack const * const kaon, StPicoTrack const * const pion);
   short const PionCharge() const ;
   short const KaonCharge() const;
   StLorentzVectorF const& lorentzVector() const;
   float const pt() const
   {
        return mD0pt;
   }
   ~StMixedD0()
   {
   }
private:
   StLorentzVectorF mLorentzVector;
   short mPionCharge;
   short mKaonCharge;
   float mD0pt;
   ClassDef(StMixedD0, 1);
};
inline short const StMixedD0::KaonCharge() const
{
    return mKaonCharge;
}
inline short const StMixedD0::PionCharge() const
{
    return mPionCharge;
}
inline StLorentzVectorF const& StMixedD0::lorentzVector() const
{
    return mLorentzVector;
}
#endif
