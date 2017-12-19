#include <limits>
#include <cmath>

#ifdef __ROOT__
#include "StD0Pion.h"
#include "StKaonPion.h"
#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"

ClassImp(StD0Pion)

//contruction function

StD0Pion::StD0Pion() : mLorentzVector{},
  mdeltaM(std::numeric_limits<float>::min()),
mD0kIdx(std::numeric_limits<unsigned short>::max()),
mD0pIdx(std::numeric_limits<unsigned short>::max()),
msPIdx(std::numeric_limits<unsigned short>::max()),
mD0CosThetaStar(std::numeric_limits<float>::max()),
mD0kaonDca(std::numeric_limits<float>::min()),
mD0pionDca(std::numeric_limits<float>::min()),
msPionDca(std::numeric_limits<float>::max()),
mD0DecayLength(std::numeric_limits<float>::min()),
mD0PointingAngle(std::numeric_limits<float>::max()) 
{}

StD0Pion::StD0Pion(StKaonPion const* aD0, StPicoTrack const& D0kaon, StPicoTrack const& D0pion, StPicoTrack const& sPion, unsigned short D0kIdx, unsigned short D0pIdx, unsigned short sPIdx, StThreeVectorF const& vtx, float bField) : StD0Pion()
{ 
    if (D0kaon.id() == sPion.id() ||
        D0pion.id() == sPion.id() )
        {
            return;
        }
    
    mD0kIdx = D0kIdx;
    mD0pIdx = D0pIdx;
    msPIdx = sPIdx;

    mD0CosThetaStar = aD0->cosThetaStar();
    mD0kaonDca = aD0->kaonDca();
    mD0pionDca = aD0->pionDca();
    mD0DecayLength = aD0->decayLength();
    mD0PointingAngle = aD0->pointingAngle();

    //softpion
    StPhysicalHelixD spHelix = sPion.dcaGeometry().helix();
    spHelix.moveOrigin(spHelix.pathLength(vtx));
    StThreeVectorF const spMom = spHelix.momentum(bField * kilogauss);
    StLorentzVectorF const spFourMom(spMom,spMom.massHypothesis(M_PION_PLUS));

    mLorentzVector = spFourMom + aD0->lorentzVector();
    mdeltaM = mLorentzVector.m()-(aD0->lorentzVector()).m();

    msPionDca = (spHelix.origin() - vtx).mag();
}

#endif

