#include <limits>
#include <cmath>

#ifdef __ROOT__
#include "StMixedD0Pion.h"
#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StMixedD0.h"
#include "StMixedSoftPion.h"

ClassImp(StMixedD0Pion)

//contruction function

StMixedD0Pion::StMixedD0Pion() : mdeltaM(std::numeric_limits<float>::min()),mLorentzVector{},mRightsign(false),mD0toPion_pt(std::numeric_limits<float>::min())
{}

StMixedD0Pion::StMixedD0Pion(StMixedD0 const& aD0, StMixedSoftPion const & aPion) : StMixedD0Pion()
{    
    mLorentzVector = aPion.lorentzVector() + aD0.lorentzVector();
    mdeltaM = mLorentzVector.m()-(aD0.lorentzVector()).m();
    mRightsign = aD0.PionCharge()*aPion.charge()>0;
    mD0toPion_pt=aD0.pt()/aPion.pt();
}


#endif