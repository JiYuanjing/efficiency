#include <limits>
#include <algorithm>

#ifdef __ROOT__
#include "StPicoKPiX.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"

ClassImp(StPicoKPiX)

StPicoKPiX::StPicoKPiX(): mKaonMomAtDca{}, mPionMomAtDca{}, mXaonMomAtDca{},
                          mKaonPionDca{std::numeric_limits<float>::max()},
                          mKaonXaonDca{std::numeric_limits<float>::max()},
                          mPionXaonDca{std::numeric_limits<float>::max()},
                          mKaonDca(std::numeric_limits<float>::min()), 
                          mPionDca(std::numeric_limits<float>::min()),
                          mXaonDca(std::numeric_limits<float>::min()),
                          mPointingAngle(std::numeric_limits<float>::max()), 
                          mDecayLength(std::numeric_limits<float>::min()),
                          mKaonIdx(std::numeric_limits<unsigned short>::max()), 
                          mPionIdx(std::numeric_limits<unsigned short>::max()),
                          mXaonIdx(std::numeric_limits<unsigned short>::max()),
                          mXaonPid(std::numeric_limits<unsigned char>::max())
{
}

//------------------------------------
StPicoKPiX::StPicoKPiX(StPicoTrack const& kaon, StPicoTrack const& pion, StPicoTrack const& xaon,
                       unsigned short const kIdx, unsigned short const pIdx, unsigned short xIdx,
                       StThreeVectorF const& vtx, float const bField, int const xaonPid) : StPicoKPiX()
{
   if (kaon.id() == pion.id() || 
       kaon.id() == xaon.id() ||
       pion.id() == xaon.id())
   {
      return;
   }

   mKaonIdx = kIdx;
   mPionIdx = pIdx;
   mXaonIdx = xIdx;
   mXaonPid = xaonPid;

   /// local variables prefixes:
   ///   k is for kaon
   ///   p is for pion
   ///   x is for xaon

   StPhysicalHelixD kHelix = kaon.dcaGeometry().helix();
   StPhysicalHelixD pHelix = pion.dcaGeometry().helix();
   StPhysicalHelixD xHelix = xaon.dcaGeometry().helix();

   // move origins of helices to the primary vertex origin
   kHelix.moveOrigin(kHelix.pathLength(vtx));
   pHelix.moveOrigin(pHelix.pathLength(vtx));
   xHelix.moveOrigin(xHelix.pathLength(vtx));

   // use straight lines approximation to get point of DCA
   StPhysicalHelixD const kStraightLine(kHelix.momentum(bField * kilogauss), kHelix.origin(), 0, kaon.charge());
   StPhysicalHelixD const pStraightLine(pHelix.momentum(bField * kilogauss), pHelix.origin(), 0, pion.charge());
   StPhysicalHelixD const xStraightLine(xHelix.momentum(bField * kilogauss), xHelix.origin(), 0, xaon.charge());

   pair<double, double> const sskp = kStraightLine.pathLengths(pStraightLine);
   pair<double, double> const sskx = kStraightLine.pathLengths(xStraightLine);
   pair<double, double> const sspx = pStraightLine.pathLengths(xStraightLine);

   StThreeVectorF const kAtDcaToP = kStraightLine.at(sskp.first);
   StThreeVectorF const pAtDcaToK = pStraightLine.at(sskp.second);
   StThreeVectorF const kAtDcaToX = kStraightLine.at(sskx.first);
   StThreeVectorF const xAtDcaToK = xStraightLine.at(sskx.second);
   StThreeVectorF const pAtDcaToX = pStraightLine.at(sspx.first);
   StThreeVectorF const xAtDcaToP = xStraightLine.at(sspx.second);

   StThreeVectorF const v0 = ( kAtDcaToP + pAtDcaToK + kAtDcaToX + xAtDcaToK + pAtDcaToX + xAtDcaToP ) / 6.;
   mKaonMomAtDca  = kHelix.momentumAt(kHelix.pathLength(v0), bField * kilogauss);
   mPionMomAtDca  = pHelix.momentumAt(pHelix.pathLength(v0), bField * kilogauss);
   mXaonMomAtDca  = xHelix.momentumAt(xHelix.pathLength(v0), bField * kilogauss);

   mKaonPionDca = (kAtDcaToP - pAtDcaToK).mag();
   mKaonXaonDca = (kAtDcaToX - xAtDcaToK).mag();
   mPionXaonDca = (pAtDcaToX - xAtDcaToP).mag();

   // calculate pointing angle and decay length
   StThreeVectorF const vtxToV0 = v0 - vtx;
   mPointingAngle = vtxToV0.angle(mKaonMomAtDca + mPionMomAtDca + mXaonMomAtDca);
   mDecayLength = vtxToV0.mag();

   // calculate DCA of tracks to primary vertex
   mKaonDca = (kHelix.origin() - vtx).mag();
   mPionDca = (pHelix.origin() - vtx).mag();
   mXaonDca = (xHelix.origin() - vtx).mag();
}

StLorentzVectorF StPicoKPiX::fourMom(double const xMassHypothesis) const
{
  return {mKaonMomAtDca+mPionMomAtDca+mXaonMomAtDca,
          mKaonMomAtDca.massHypothesis(M_KAON_MINUS)+
          mPionMomAtDca.massHypothesis(M_PION_PLUS)+
          mXaonMomAtDca.massHypothesis(xMassHypothesis)};
}

StLorentzVectorF StPicoKPiX::kaonPionFourMom() const
{
  return {mKaonMomAtDca+mPionMomAtDca, 
          mKaonMomAtDca.massHypothesis(M_KAON_MINUS)+
          mPionMomAtDca.massHypothesis(M_PION_PLUS)};
}

StLorentzVectorF StPicoKPiX::kaonXaonFourMom(double const xMassHypothesis) const
{
  return {mKaonMomAtDca+mXaonMomAtDca, 
          mKaonMomAtDca.massHypothesis(M_KAON_MINUS)+
          mXaonMomAtDca.massHypothesis(xMassHypothesis)};
}

StLorentzVectorF StPicoKPiX::pionXaonFourMom(double const xMassHypothesis) const
{
  return {mPionMomAtDca+mXaonMomAtDca, 
          mPionMomAtDca.massHypothesis(M_KAON_MINUS)+
          mXaonMomAtDca.massHypothesis(xMassHypothesis)};
}

float StPicoKPiX::dcaDaughters() const
{
  return std::max({mKaonPionDca,mKaonXaonDca,mPionXaonDca});
}

#endif // __ROOT__
