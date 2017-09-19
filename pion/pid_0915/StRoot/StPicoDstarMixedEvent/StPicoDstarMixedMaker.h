#ifndef StPicoDstarMixedMaker_h
#define StPicoDstarMixedMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "StDstarMixedEvent.h"
class TString;
class TFile;
class TNtuple;
class StPicoEvent;
class StPicoD0Event;
class StKaonPion;
class StD0Pion;
class StPicoTrack;
class StPicoDstMaker;
class StPicoDstarMixedHists;
class StRefMultCorr;

class StPicoDstarMixedMaker : public StMaker
{
public:
   StPicoDstarMixedMaker(char const * name, TString const inputFilesList,
                    TString const outBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil);
    virtual ~StPicoDstarMixedMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;
    void fillQaHistograms(bool b = true);
    void fillBackgroundTrees(bool b = true);
    void addmixedevent(bool b=true)
    {
        mMixedEvent = b;
    }
    void fillSoftPionEff(bool b = true )
    {
        mFillSoftPionEff = b;
    }
    void reconstructD(bool b = true)
    {
        mReconstructD = b;
    }
    void softpionQa(bool b = true)
    {
        mSoftPionQa = b;
    }
private:

    StPicoDstarMixedMaker() {}
    void readNextEvent();

    int  getD0PtIndex(StKaonPion const* ) const;
    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodEvent(StPicoEvent const*, StThreeVectorF const& vtx) const;
    bool isGoodQaTrack(StPicoTrack const* ,StThreeVectorF const& momentum ,double dca) const;
    bool isGoodTrack(StPicoTrack const*, StThreeVectorF const&) const;
    bool isGoodSoftPionTrack(StPicoTrack const* const trk, StThreeVectorF const& mom, float spDca) const;
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isTpcProton(StPicoTrack const*) const;
    bool isTofPion(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofKaon(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofProton(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isGoodPair(StKaonPion const*) const;
    bool isD0SideBand(float m) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const& vtx) const;
    int trkHalf(StPicoTrack const*, StThreeVectorF const& vtx) const;
    bool isGoodD0(StKaonPion const*) const;
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StRefMultCorr* mGRefMultCorrUtil;
 
    TString mInputFilesList;
    TString mOutFileBaseName;
    TChain* mChain;
    int mEventCounter;
    bool mFillQaHists;
    bool mFillBackgroundTrees;  
    bool mFillSoftPionEff;
    bool mReconstructD;
    bool mSoftPionQa;
   // -------------- USER variables -------------------------
   // add your member variables here. 
   // Remember that ntuples size can be really big, use histograms where appropriate 

   StPicoDstarMixedHists* mHists;

//mixed event varibles
public:
void setEventsBufferSize(int bufferSize);
//void mixEvents();

private:
    bool mMixedEvent;
    int mEventsBufferSize;
    std::vector <StDstarMixedEvent> mEvents[10][9];

    int mCentBin, mVzBin;
   ClassDef(StPicoDstarMixedMaker, 1)
};

inline void StPicoDstarMixedMaker::fillQaHistograms(bool b) { mFillQaHists = b;}
inline void StPicoDstarMixedMaker::fillBackgroundTrees(bool b) { mFillBackgroundTrees = b;}
inline int StPicoDstarMixedMaker::getEntries() const
{
   return mChain ? mChain->GetEntries() : 0;
}
inline void StPicoDstarMixedMaker::readNextEvent()
{
   mChain->GetEntry(mEventCounter++);
}
//used to set the mixed event buffersize//
inline void StPicoDstarMixedMaker::setEventsBufferSize(int bufferSize)
{
   mEventsBufferSize = bufferSize;
}
#endif
