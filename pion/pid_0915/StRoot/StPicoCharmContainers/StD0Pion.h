#ifndef StD0Pion_hh
#define StD0Pion_hh
#ifdef __ROOT__

#include <cmath>

#include "TObject.h"
#include "TClonesArray.h"
#include "StLorentzVectorF.hh"

class StPicoTrack;
class StPicoEvent;
class StKaonPion;
class StD0Pion : public TObject
{
    public:
    StD0Pion();
    StD0Pion(StKaonPion const* aD0, StPicoTrack const& D0kaon, StPicoTrack const& D0pion, StPicoTrack const& sPion, unsigned short D0kIdx, unsigned short D0pIdx, unsigned short sPIdx, StThreeVectorF const& vtx, float bField);
    ~StD0Pion() {}

    StLorentzVectorF const & lorentzVector() const { return mLorentzVector;}
//D_star
    float deltaM() const { return mdeltaM;}
    float m() const { return mLorentzVector.m();}
    float pt() const { return mLorentzVector.perp();}
    float eta() const { return mLorentzVector.pseudoRapidity();}
    float phi() const { return mLorentzVector.phi();}
//    StLorentzVector fourMom() const;

    float D0kaonDca() const{ return mD0kaonDca;}
    float D0pionDca() const{ return mD0pionDca;}
    float sPionDca() const{ return msPionDca;}

    unsigned short   D0kIdx() const{ return  mD0kIdx;}
    unsigned short   D0pIdx() const{ return mD0pIdx;}
    unsigned short   sPIdx() const{ return msPIdx;}
    unsigned short D0DecayLength() const{ return mD0DecayLength;}
    unsigned short D0PointingAngle() const{ return mD0PointingAngle;}

    float D0CosThetaStar() const { return mD0CosThetaStar;}
    float D0perpDcaToVtx() const { return mD0DecayLength*std::sin(mD0PointingAngle);}

private:
    float mdeltaM;
    StLorentzVectorF mLorentzVector;
    unsigned short mD0kIdx;
    unsigned short mD0pIdx;
    unsigned short msPIdx;
    
    float mD0CosThetaStar;
    float mD0kaonDca;
    float mD0pionDca;
    float msPionDca;
    float mD0DecayLength;
    float mD0PointingAngle;    

    ClassDef(StD0Pion,1)
};
#endif
#endif
