#include "StMixedSoftPion.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StThreeVectorF.hh"
#include <limits>
#include "StLorentzVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
ClassImp(StMixedSoftPion)
StMixedSoftPion::StMixedSoftPion() : mOrigin(StThreeVectorF()), mMom(StThreeVectorF()), mCharge(std::numeric_limits<short>::min())
{
}
StMixedSoftPion::StMixedSoftPion(StThreeVectorF const & pVtx, float B, StPicoTrack const* const picoTrack)
{
   StPhysicalHelixD helix = picoTrack->helix();
   helix.moveOrigin(helix.pathLength(pVtx));
   mOrigin = helix.origin();
   mMom = helix.momentum(B * kilogauss);
   mCharge = picoTrack->charge();
   StLorentzVectorF pionfourmom(mMom,mMom.massHypothesis(M_PION_PLUS));
   mLorentzVector = pionfourmom;
   mPt = mLorentzVector.perp();
  //mPt = mMom.perp();
}
StMixedSoftPion::StMixedSoftPion(StMixedSoftPion const * t) : mOrigin(t->mOrigin), mMom(t->mMom), mCharge(t->mCharge), mLorentzVector(t->mLorentzVector)
{
}

