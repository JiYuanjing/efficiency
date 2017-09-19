#include "StMixedD0.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StThreeVectorF.hh"
#include <limits>
#include "StLorentzVectorF.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
ClassImp(StMixedD0)
StMixedD0::StMixedD0() : mPionCharge(std::numeric_limits<short>::min()),mKaonCharge(std::numeric_limits<short>::min()), mLorentzVector(StLorentzVectorF()),mD0pt(std::numeric_limits<float>::min())
{
}
StMixedD0::StMixedD0(StThreeVectorF const & pVtx, float B, StKaonPion const * const aD0,StPicoTrack const* const kaon, StPicoTrack const* const pion)
{
   mLorentzVector = aD0->lorentzVector();
   mPionCharge = pion->charge();
   mKaonCharge = kaon->charge();
   mD0pt = mLorentzVector.perp();
}
StMixedD0::StMixedD0(StMixedD0 const * t) : mLorentzVector(t->mLorentzVector), mPionCharge(t->mPionCharge), mKaonCharge(t->mKaonCharge),mD0pt(t->mD0pt)
{
}
