#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"

#include "StPicoKPiXEvent.h"
#include "StPicoKPiX.h"

ClassImp(StPicoKPiXEvent)

TClonesArray *StPicoKPiXEvent::fgKaonPionXaonArray = nullptr;


StPicoKPiXEvent::StPicoKPiXEvent() : mRunId(-1), mEventId(-1), mNKaonPionXaon(0), mKaonPionXaonArray(nullptr)
{
   if (!fgKaonPionXaonArray) fgKaonPionXaonArray = new TClonesArray("StPicoKPiX");
   mKaonPionXaonArray = fgKaonPionXaonArray;
}


void StPicoKPiXEvent::addPicoEvent(StPicoEvent const& picoEvent)
{
   // StPicoEvent variables
   mRunId = picoEvent.runId();
   mEventId = picoEvent.eventId();
}


void StPicoKPiXEvent::clear(char const *option)
{
   mKaonPionXaonArray->Clear(option);
   mRunId = -1;
   mEventId = -1;
   mNKaonPionXaon = 0;
}

void StPicoKPiXEvent::addKPiX(StPicoKPiX const& t)
{
   new((*mKaonPionXaonArray)[mNKaonPionXaon++]) StPicoKPiX(t);
}
