#ifndef StMixedD0Pion_hh
#define StMixedD0Pion_hh
#ifdef __ROOT__

#include <cmath>

#include "TObject.h"
#include "TClonesArray.h"
#include "StLorentzVectorF.hh"
class StMixedD0;
class StMixedSoftPion;
class StMixedD0Pion : public TObject
{
    public:
    StMixedD0Pion();
    StMixedD0Pion(StMixedD0 const& aD0, StMixedSoftPion const  & aPion);
    ~StMixedD0Pion() {}

    StLorentzVectorF const & lorentzVector() const { return mLorentzVector;}
//D_star
    float deltaM() const { return mdeltaM;}
    float m() const { return mLorentzVector.m();}
    float pt() const { return mLorentzVector.perp();}
    float eta() const { return mLorentzVector.pseudoRapidity();}
    float phi() const { return mLorentzVector.phi();}
    float D0toPion_pt() const{
        return mD0toPion_pt;
    }
    bool rightsign() const { return mRightsign;};

private:
    float mdeltaM;
    StLorentzVectorF mLorentzVector;
    bool mRightsign;
    float mD0toPion_pt;
    ClassDef(StMixedD0Pion,1)
};
#endif
#endif
