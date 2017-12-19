#include <limits>
#include <cmath>

#ifdef __ROOT__
#include "StKaonPion.h"

#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"

ClassImp(StKaonPion)

//initialize construction function
StKaonPion::StKaonPion(): mLorentzVector{},
                          mPointingAngle(std::numeric_limits<float>::max()), 
                          mDecayLength(std::numeric_limits<float>::min()),
                          mKaonDca(std::numeric_limits<float>::min()), 
                          mPionDca(std::numeric_limits<float>::min()),
                          mKaonIdx(std::numeric_limits<unsigned short>::max()), 
                          mPionIdx(std::numeric_limits<unsigned short>::max()),
                          mDcaDaughters(std::numeric_limits<float>::max()), 
                          mCosThetaStar(std::numeric_limits<float>::max())
{
}

StKaonPion::StKaonPion(StPicoTrack const& kaon, StPicoTrack const& pion,
                       unsigned short const kIdx, unsigned short const pIdx,
                       StThreeVectorF const& vtx, float const bField) : StKaonPion()
{
   if (kaon.id() == pion.id())
   {
      mKaonIdx = std::numeric_limits<unsigned short>::quiet_NaN();
      mPionIdx = std::numeric_limits<unsigned short>::quiet_NaN();
      return;
   }

   mKaonIdx = kIdx;
   mPionIdx = pIdx;

   /// prefixes code:
   ///   k means kaon
   ///   p means pion
   ///   kp means kaon-pion pair


   // to be used for testing with preview II pico production
   StPhysicalHelixD kHelix = kaon.dcaGeometry().helix();
   StPhysicalHelixD pHelix = pion.dcaGeometry().helix();

   // move origins of helices to the primary vertex origin
   kHelix.moveOrigin(kHelix.pathLength(vtx));
   pHelix.moveOrigin(pHelix.pathLength(vtx));

   // use straight lines approximation to get point of DCA of kaon-pion pair
   StThreeVectorF const kMom = kHelix.momentum(bField * kilogauss);
   StThreeVectorF const pMom = pHelix.momentum(bField * kilogauss);
   StPhysicalHelixD const kStraightLine(kMom, kHelix.origin(), 0, kaon.charge());
   StPhysicalHelixD const pStraightLine(pMom, pHelix.origin(), 0, pion.charge());

   pair<double, double> const ss = kStraightLine.pathLengths(pStraightLine);
   StThreeVectorF const kAtDcaToPion = kStraightLine.at(ss.first);
   StThreeVectorF const pAtDcaToKaon = pStraightLine.at(ss.second);

   // calculate DCA of pion to kaon at their DCA
   mDcaDaughters = (kAtDcaToPion - pAtDcaToKaon).mag();

   // calculate Lorentz vector of kaon-pion pair
   StThreeVectorF const kMomAtDca = kHelix.momentumAt(ss.first, bField * kilogauss);
   StThreeVectorF const pMomAtDca = pHelix.momentumAt(ss.second, bField * kilogauss);

   StLorentzVectorF const kFourMom(kMomAtDca, kMomAtDca.massHypothesis(M_KAON_PLUS));
   StLorentzVectorF const pFourMom(pMomAtDca, pMomAtDca.massHypothesis(M_PION_PLUS));

   mLorentzVector = kFourMom + pFourMom;

   // calculate cosThetaStar
   StLorentzVectorF const kpFourMomReverse(-mLorentzVector.px(), -mLorentzVector.py(), -mLorentzVector.pz(), mLorentzVector.e());
   StLorentzVectorF const kFourMomStar = kFourMom.boost(kpFourMomReverse);
   mCosThetaStar = std::cos(kFourMomStar.vect().angle(mLorentzVector.vect()));

   // calculate pointing angle and decay length
   StThreeVectorF const vtxToV0 = (kAtDcaToPion + pAtDcaToKaon) * 0.5 - vtx;
   mPointingAngle = vtxToV0.angle(mLorentzVector.vect());
   mDecayLength = vtxToV0.mag();

   // calculate DCA of tracks to primary vertex
   mKaonDca = (kHelix.origin() - vtx).mag();
   mPionDca = (pHelix.origin() - vtx).mag();
}
#endif // __ROOT__
