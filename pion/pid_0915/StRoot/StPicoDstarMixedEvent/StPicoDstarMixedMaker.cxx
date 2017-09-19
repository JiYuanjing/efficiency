 /* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoCharmContainers/StPicoD0Event.h"
#include "StPicoCharmContainers/StKaonPion.h"
#include "StPicoCharmContainers/StD0Pion.h"
#include "StPicoDstarMixedMaker.h"
#include "StPicoDstarMixedHists.h"
#include "StMixedD0.h"
#include "StMixedSoftPion.h"
#include "StDstarMixedEvent.h"
#include "StMixedD0Pion.h"
#include "StAnaCuts.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StMemStat.h"

ClassImp(StPicoDstarMixedMaker)

StPicoDstarMixedMaker::StPicoDstarMixedMaker(char const * name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil):
   StMaker(name), mPicoDstMaker(picoDstMaker), mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
   mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),
   mChain(NULL), mEventCounter(0), mFillQaHists(false), mFillBackgroundTrees(false), mHists(NULL),mMixedEvent(false),mFillSoftPionEff(true),mReconstructD(true), mEventsBufferSize(6)
{}

Int_t StPicoDstarMixedMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFilesList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoDstarMixedMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoDstarMixedMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
   // -------------- USER VARIABLES -------------------------
   mHists = new StPicoDstarMixedHists(mOutFileBaseName, mReconstructD,mFillQaHists, mFillBackgroundTrees, mFillSoftPionEff, mSoftPionQa);
   setEventsBufferSize(6);
   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarMixedMaker::~StPicoDstarMixedMaker()
{

}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Finish()
{
   mHists->closeFile();
   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Make()
{
   StMemStat mem;
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoDstarMixedMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoDstarMixedMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if (mPicoD0Event->runId() != picoDst->event()->runId() ||
         mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
      LOG_ERROR << " StPicoDstarMixedMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << "\n";
      LOG_ERROR << " StPicoDstarMixedMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync." << endm;
      exit(1);
   }

   // -------------- USER ANALYSIS -------------------------

   mGRefMultCorrUtil->init(picoDst->event()->runId());

   if (!mGRefMultCorrUtil)
   {
      LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
      return kStWarn;
   }

   if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId()))
   {
      //cout<<"This is a bad run from mGRefMultCorrUtil! Skip! " << endl;
      return kStOK;
   }

   if (isGoodTrigger(picoDst->event()))
   {
      mHists->addEventBeforeCut(picoDst->event());

      StThreeVectorF pVtx = picoDst->event()->primaryVertex();
      if (isGoodEvent(picoDst->event(), pVtx))
      {
         TClonesArray const* aKaonPion = mPicoD0Event->kaonPionArray();
         mHists->addEvent(picoDst->event());

         mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(), pVtx.z(), picoDst->event()->ZDCx()) ;

         int centrality  = mGRefMultCorrUtil->getCentralityBin9();
         const double reweight = mGRefMultCorrUtil->getWeight();
         const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();
         mHists->addCent(refmultCor, centrality, reweight, pVtx.z());
         //Basiclly add some QA plots
         UInt_t nTracks = picoDst->numberOfTracks();
          //vz range is (-6,6)
          int const vz_bin = (int)((6 + pVtx.z()) / 1.2);
         if (centrality < 0 || centrality > 8) return kStOk;
         if (vz_bin < 0  ||  vz_bin > 9) return kStOk;

         //help to efficiency correction, used for D0//
         if (mFillQaHists)
         {
            for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
            {
               StPicoTrack const* trk = picoDst->track(iTrack);
               if (!trk) continue;
               StPhysicalHelixD helix = trk->helix();
               float dca = float(helix.geometricSignedDistance(pVtx));
               StThreeVectorF momentum = trk->gMom(pVtx, picoDst->event()->bField());

               if (!isGoodQaTrack(trk, momentum, dca)) continue;

               StThreeVectorF dcaPoint = helix.at(helix.pathLength(pVtx.x(), pVtx.y()));
               float dcaZ = dcaPoint.z() - pVtx.z();
               double dcaXy = helix.geometricSignedDistance(pVtx.x(), pVtx.y());

               bool tpcPion = isTpcPion(trk);
               bool tpcKaon = isTpcKaon(trk);
               bool tpcProton = isTpcProton(trk);
               float piBeta = getTofBeta(trk, pVtx);
               float kBeta = piBeta;
               float pBeta = piBeta;
               bool piTofAvailable = !isnan(piBeta) && piBeta > 0;
               bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
               bool pTofAvailable = !isnan(pBeta) && pBeta > 0;
               bool tofPion = isTofPion(trk, piBeta, pVtx);
               bool tofKaon = isTofKaon(trk, kBeta, pVtx);
               bool tofProton = isTofProton(trk, pBeta, pVtx);

               bool goodPion = (piTofAvailable && tofPion && tpcPion) || (!piTofAvailable && tpcPion);//Always require TPC
               bool goodKaon = (kTofAvailable && tofKaon && tpcKaon) || (!kTofAvailable && tpcKaon);
               bool goodProton = (pTofAvailable && tofProton && tpcProton) || (!pTofAvailable && tpcProton);
               // bool goodKaon = (momentum.perp() <= 1.6 && kTofAvailable && tofKaon && tpcKaon) || (momentum.perp() > 1.6 && tpcKaon);//strict Kaon pid

               // cout<<"mem.Used() = " << mem.Used()<<endl;
               // cout<<"mem.Free() = " << mem.Free()<<endl;
               // cout<<"mem.ProgSize() = " << mem.ProgSize()<<endl;

               // bool goodPion = piTofAvailable && tofPion && tpcPion;//strict Pion pid
               // bool goodKaon = kTofAvailable && tofKaon && tpcKaon;//strict Kaon pid
               // bool goodProton = pTofAvailable && tofProton && tpcProton;//strict Proton pid

               if (trk  && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton))
               {
                  mHists->addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //add Dca distribution
               }
               if (trk  && fabs(dca) < 1.5 && (goodPion || goodKaon || goodProton))
               {
                  mHists->addTpcDenom1(goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add Tpc Denominator
               }
               if (trk && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton) && fabs(dcaXy) < 1. && fabs(dcaZ) < 1.)
               {
                  mHists->addHFTNumer1(goodPion, goodKaon, goodProton, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add HFT Numerator
               }
            } // .. end tracks loop
         }

//mixed event
          StDstarMixedEvent event(pVtx, picoDst->event()->bField(), reweight);
         if (mReconstructD)
         {
         for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
         {
            StKaonPion const* kp = (StKaonPion*)aKaonPion->UncheckedAt(idx);

            bool goodPair = isGoodPair(kp);
            if (!goodPair && !mFillBackgroundTrees) continue;
            StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
            StPicoTrack const* pion = picoDst->track(kp->pionIdx());

            if (!isGoodTrack(kaon, pVtx) || !isGoodTrack(pion, pVtx)) continue;
  
            // PID
            if (!isTpcPion(pion) || !isTpcKaon(kaon)) continue;
            float piBeta = getTofBeta(pion, pVtx);
            float kBeta = getTofBeta(kaon, pVtx);
            bool piTofAvailable = !isnan(piBeta) && piBeta > 0;
            bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
            bool tofPion = piTofAvailable ? isTofPion(pion, piBeta, pVtx) : true;//this is bybrid pid, not always require tof
            bool tofKaon = kTofAvailable ? isTofKaon(kaon, kBeta, pVtx) : true;//this is bybrid pid, not always require tof

            bool tof = tofPion && tofKaon;

            bool D0unlike = kaon->charge() * pion->charge() < 0 ? true : false;
//D0
            if (goodPair) mHists->addKaonPion(kp, D0unlike, true, tof, centrality, reweight);

            if (mFillBackgroundTrees && tof)
            {
               if (centrality < 0) continue;
               int const ptBin = getD0PtIndex(kp);

               if (D0unlike)
               {
                  if (kp->m() > anaCuts::likeSignMassRange.first && kp->m() < anaCuts::likeSignMassRange.second) mHists->addBackground(kp, kaon, pion, ptBin, false);
               }
               else if (isD0SideBand(kp->m()))
               {
                  mHists->addBackground(kp, kaon, pion, ptBin, true);
               }
            }
     //Dstar 
            if (!tof) continue;
            if (!goodPair) continue;
            if (!D0unlike) continue;

//several judgement condition//
//PID(not require tof match), idx^_^, charge(K+pi-pi- or K-pi+pi+), the pt of softpion, dca, pt(D0)/pt(Dstar) //
//here is left for sideband background
            if (!( isGoodD0(kp) || isD0SideBand(kp->m()) )) continue;
//add softpion loop to reconstruct Dstar//
//mixed event, add D0
            if (mReconstructD && mMixedEvent && isGoodD0(kp)){
                StMixedD0 mixD0(pVtx,picoDst->event()->bField(), kp, kaon, pion);
                event.addD0(mixD0);
              }
//pair the D0 with a softpion
//this loop is inside the D0 loop, so need to loop it again seperately outside the D0 loop if mixed event is required.
            for (unsigned short iTrack = 0; iTrack<nTracks; ++iTrack)
            {
               StPicoTrack* trk = picoDst->track(iTrack);
               if (!trk) continue;
               if (iTrack==kp->kaonIdx() ||
                   iTrack==kp->pionIdx())  continue;

               StThreeVectorF spMom = trk->gMom(pVtx, picoDst->event()->bField());
               float spDca = float((trk->helix()).geometricSignedDistance(pVtx));
               if (!isGoodSoftPionTrack(trk, spMom, spDca)) continue;
//PID 
               if (!isTpcPion(trk)) continue;
               float spiBeta = getTofBeta(trk, pVtx);
               bool spiTofAvailable = !isnan(spiBeta) && spiBeta > 0;
               bool tofSoftPion = spiTofAvailable ? isTofPion(trk, spiBeta, pVtx) : true;//this is bybrid pid, not always require tof

               if (!tofSoftPion) continue;

              // if (kp->pt()/spMom.perp()>20 || kp->pt()/spMom.perp()<7) continue;
               StD0Pion D0Pion(kp,*kaon , *pion, *trk, kp->kaonIdx(),kp->pionIdx(),iTrack,pVtx, picoDst->event()->bField());
               bool rightsign = pion->charge() * trk->charge() >0 ? true : false;
               if ((D0Pion.m()-kp->m())<0.1 || (D0Pion.m()-kp->m())>0.2) continue;  
               if (isGoodD0(kp)) { 
                 mHists->addD0SoftPion(&D0Pion,kp,rightsign,centrality,reweight);
                 if (mSoftPionQa){
                  if (D0Pion.deltaM()<0.144 && D0Pion.deltaM()<0.147)
                   mHists->addSoftPionQa(spMom,spDca, centrality, reweight,trk->charge());
                 }
              }
               if (isD0SideBand(kp->m()))
               mHists->addSideBandBackground(&D0Pion,kp,rightsign,centrality,reweight);
           }//end of iTrack loop
         } // end of kaonPion loop
        }//reconstructD 

//caculate the softpion eff
         if (mFillSoftPionEff)
         {
          for (unsigned short iTrack = 0; iTrack<nTracks; ++iTrack)
          {
             StPicoTrack* trk = picoDst->track(iTrack);
             if (!trk) continue;

             StThreeVectorF spMom = trk->gMom(pVtx, picoDst->event()->bField());
             float spDca = float((trk->helix()).geometricSignedDistance(pVtx));
             if (!isGoodSoftPionTrack(trk, spMom, spDca)) continue;
             
             double diffInvBeta = -999;
             double nSigmaPion = -999;
             double nSigmaKaon = -999;
             double nSigmaProton = -999;
            //softpion PID 
           
            nSigmaPion = trk->nSigmaPion();
            nSigmaKaon = trk->nSigmaKaon();
            nSigmaProton = trk->nSigmaProton();
            bool isTpcPion = fabs(nSigmaPion) < anaCuts::nSigmaPion && fabs(nSigmaKaon)>4 && fabs(nSigmaProton)>4;
            float spiBeta = getTofBeta(trk, pVtx);
            bool spiTofAvailable = !isnan(spiBeta) && spiBeta > 0;
            
            bool tofPion = false;
            bool nottofkaon = false;
            bool nottofproton = false;
            float beta_pi = -999;
            float beta_k = -999;
            float beta_pr = -999;
            if (spiTofAvailable){
                 double ptot = trk->gMom(pVtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
                 beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
                 beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
                 beta_pr = ptot / sqrt(ptot * ptot + M_PROTON * M_PROTON);
                 diffInvBeta = 1.0 / spiBeta - 1.0 / beta_pi;
                 nottofkaon = fabs(1.0 / spiBeta - 1.0 / beta_k)>0.04;
                 nottofproton = fabs(1.0 / spiBeta - 1.0 / beta_pr)>0.04;
                 tofPion = fabs(diffInvBeta) < anaCuts::piTofBetaDiff && nottofproton && nottofkaon;
            }
            mHists->addSoftPionEff(spMom,isTpcPion,spiTofAvailable, tofPion, diffInvBeta, nSigmaPion, beta_pi, centrality, reweight, trk->charge());           
          }//itrack
         }//softpioneff

//mixed event, loop and add softpion.
if (mReconstructD&&mMixedEvent){
  //if no D0, then delete the event.
            if (event.getNoD0s()<1) {
              return kStOK;
            } 

           for (unsigned short iTrack = 0; iTrack<nTracks; ++iTrack)
          {
               StPicoTrack* trk = picoDst->track(iTrack);
               if (!trk) continue;
               StThreeVectorF spMom = trk->gMom(pVtx, picoDst->event()->bField());
               float spDca = float((trk->helix()).geometricSignedDistance(pVtx));
               if (!isGoodSoftPionTrack(trk, spMom, spDca)) continue;
//PID 
               if (!isTpcPion(trk)) continue;
               float spiBeta = getTofBeta(trk, pVtx);
               bool spiTofAvailable = !isnan(spiBeta) && spiBeta > 0;
               bool tofSoftPion = spiTofAvailable ? isTofPion(trk, spiBeta, pVtx) : true;//this is bybrid pid, not always require tof
               if (!tofSoftPion) continue;

               StMixedSoftPion mixsoftpion(pVtx, picoDst->event()->bField(), trk);
               event.addSoftPion(mixsoftpion); 
            }

           if (event.getNoPions()<1) {
              return kStOK;
            }

//add new event
          mEvents[vz_bin][centrality].push_back(event);
          //if buffer is full, then mixed the event
          int mEventsize = mEvents[vz_bin][centrality].size();
          mHists->addeventsinbuffer(vz_bin,centrality,mEventsize-1);
          if (mEventsize==mEventsBufferSize)
          {
            //mixevent with each event in the buffer
            //the last event in the buffer is the current event!
            int ND0_1 = mEvents[vz_bin][centrality][mEventsize-1].getNoD0s();
            int NPion_1 = mEvents[vz_bin][centrality][mEventsize-1].getNoPions();
            int ND0_i;
            int NPion_i;
            for (int iEvt=0; iEvt<mEventsBufferSize-1;iEvt++){
              NPion_i=mEvents[vz_bin][centrality][iEvt].getNoPions();
              ND0_i=mEvents[vz_bin][centrality][iEvt].getNoD0s();

              for (int jD0_1=0;jD0_1<ND0_1;jD0_1++){
                for (int jPion_i=0;jPion_i<NPion_i;jPion_i++){
                StMixedD0Pion mixedD0Pion(mEvents[vz_bin][centrality][mEventsize-1].D0At(jD0_1),mEvents[vz_bin][centrality][iEvt].pionAt(jPion_i));
              //  if (mixedD0Pion.D0toPion_pt()< anaCuts::D0toPion_pt.first|| mixedD0Pion.D0toPion_pt()>anaCuts::D0toPion_pt.second) continue;
                mHists->addMixedEventBackground(mixedD0Pion, mixedD0Pion.rightsign(),centrality,reweight);
                }
              }
              
              for (int jPion_1=0;jPion_1<NPion_1;jPion_1++){
                for (int jD0_i=0;jD0_i<ND0_i;jD0_i++){
                StMixedD0Pion mixedD0Pion(mEvents[vz_bin][centrality][iEvt].D0At(jD0_i),mEvents[vz_bin][centrality][mEventsize-1].pionAt(jPion_1));
              //  if (mixedD0Pion.D0toPion_pt()< anaCuts::D0toPion_pt.first|| mixedD0Pion.D0toPion_pt()>anaCuts::D0toPion_pt.second) continue;
                mHists->addMixedEventBackground(mixedD0Pion, mixedD0Pion.rightsign(),centrality,reweight);
                }
              }
            }//eventsbuffer loop

          //remove the oldest event
      
            mEvents[vz_bin][centrality].erase( mEvents[vz_bin][centrality].begin());
          } //mEventsize>=mEventsBufferSize
         } //end of mixed event
        } // end of isGoodEvent
       } // end of isGoodTrigger

   return kStOK;
}
//-----------------------------------------------------------------------------
int StPicoDstarMixedMaker::getD0PtIndex(StKaonPion const* const kp) const
{
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((kp->pt() >= anaCuts::PtBinsEdge[i]) && (kp->pt() < anaCuts::PtBinsEdge[i + 1]))
         return i;
   }
   return anaCuts::nPtBins - 1;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent, StThreeVectorF const& pVtx) const
{
   return fabs(pVtx.z()) < anaCuts::vz &&
          fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
          !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror && fabs(pVtx.z()) < anaCuts::Verror) &&
          sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vrcut;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
   for (auto trg : anaCuts::triggers)
   {
      if (picoEvent->isTrigger(trg)) return true;
   }

   return false;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isGoodQaTrack(StPicoTrack const* const trk, StThreeVectorF const& momentum, double const dca) const
{
   return trk->gPt() > anaCuts::qaGPt && trk->nHitsFit() >= anaCuts::qaNHitsFit && fabs(momentum.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());


   return mom.perp() > anaCuts::minPt &&
          trk->nHitsFit() >= anaCuts::nHitsFit &&
          fabs(mom.pseudoRapidity()) <= anaCuts::Eta;
}
//----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isGoodSoftPionTrack(StPicoTrack const* const trk, StThreeVectorF const& mom, float spDca) const
{
   bool bnHistsFit = trk->nHitsFit() >= anaCuts::nHitsFit;
   bool bPt = mom.perp()>anaCuts::ptSoftPion_min;
   bool bDca =std::fabs(spDca)<=anaCuts::DcaSoftPion;
   bool bEta = std::fabs(mom.pseudoRapidity())<=anaCuts::Eta;
   return bnHistsFit && bPt && bDca && bEta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isTpcPion(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isTpcKaon(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isTpcProton(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaProton()) < anaCuts::nSigmaProton;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isGoodPair(StKaonPion const* const kp) const
{
   int tmpIndex = getD0PtIndex(kp);
   return cos(kp->pointingAngle()) > anaCuts::cosTheta[tmpIndex] &&
          kp->pionDca() > anaCuts::pDca[tmpIndex] && kp->kaonDca() > anaCuts::kDca[tmpIndex] &&
          kp->dcaDaughters() < anaCuts::dcaDaughters[tmpIndex] &&
          kp->decayLength() > anaCuts::decayLength[tmpIndex] &&
          fabs(kp->lorentzVector().rapidity()) < anaCuts::RapidityCut &&
          ((kp->decayLength()) * sin(kp->pointingAngle())) < anaCuts::dcaV0ToPv[tmpIndex];
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isD0SideBand(float const m) const
{
  bool sideband=false;
  if (m< anaCuts::sideBandMassRange0.first&& m> anaCuts::sideBandMassRange0.second) sideband=true;
  if (m> anaCuts::sideBandMassRange0.first&& m<anaCuts::sideBandMassRange0.second) sideband = true;
    return sideband;
   }
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isTofKaon(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofKaon = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
      tofKaon = fabs(1 / beta - 1 / beta_k) < anaCuts::kTofBetaDiff ? true : false;
   }

   return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isTofPion(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofPion = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < anaCuts::piTofBetaDiff ? true : false;
   }

   return tofPion;
}
//-----------------------------------------------------------------------------
bool StPicoDstarMixedMaker::isTofProton(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofProton = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_p = ptot / sqrt(ptot * ptot + M_PROTON * M_PROTON);
      tofProton = fabs(1 / beta - 1 / beta_p) < anaCuts::pTofBetaDiff ? true : false;
   }

   return tofProton;
}
//-----------------------------------------------------------------------------
float StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0)
   {
      StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if (tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            StPhysicalHelixD helix = trk->helix();
            float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
//-----------------------------------------------------------------------------
int StPicoDstarMixedMaker::trkHalf(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());
   if (mom.phi() > anaCuts:: rightHalfLowEdge && mom.phi() < anaCuts:: rightHalfHighEdge) return +1; //right side
   else return -1;//lest side

}
//left to add some cut condition, no use now//
bool StPicoDstarMixedMaker::isGoodD0(StKaonPion const* const kp) const
{
            bool mass = (kp->m()<anaCuts::mD0_max) && (kp->m()>anaCuts::mD0_min);
//            bool thetastar = (kp->cosThetaStar()<anaCuts::cosThetaStar);
//           bool pt = (kp->pt()>anaCuts::ptD0_min);
           return mass; // && thetastar && pt;
}
