#ifndef StPicoKPiXEvent__h
#define StPicoKPiXEvent__h

/* **************************************************
 *  A specialized class for storing eventwise KPiX
 *  candidates. 
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include <cstddef>

#include "TObject.h"
#include "TClonesArray.h"

class StPicoEvent;
class StPicoKPiX;

class StPicoKPiXEvent : public TObject
{
public:
   StPicoKPiXEvent();
   ~StPicoKPiXEvent(){ clear("C");}
   void clear(char const *option = "");
   void addPicoEvent(StPicoEvent const& picoEvent);
   void addKPiX(StPicoKPiX const&);

   Int_t runId()   const;
   Int_t eventId() const;
   int   nKaonPionXaon()  const;
   TClonesArray const* kaonPionXaonArray()   const;

private:
   // some variables below are kept in ROOT types to match the same ones in StPicoEvent
   Int_t mRunId;           // run number
   Int_t mEventId;         // event number
   int   mNKaonPionXaon;

   TClonesArray*        mKaonPionXaonArray;
   static TClonesArray* fgKaonPionXaonArray;

   ClassDef(StPicoKPiXEvent, 1)
};
inline TClonesArray const* StPicoKPiXEvent::kaonPionXaonArray()   const { return mKaonPionXaonArray;}
inline Int_t StPicoKPiXEvent::runId()   const { return mRunId; }
inline Int_t StPicoKPiXEvent::eventId() const { return mEventId; }
inline int   StPicoKPiXEvent::nKaonPionXaon()  const { return mNKaonPionXaon;}
#endif
