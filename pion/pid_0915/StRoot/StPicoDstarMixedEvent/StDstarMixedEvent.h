#ifndef StDstarMixedEvent_hh
#define StDstarMixedEvent_hh
#include <math.h>
#include <vector>
#include "TVector2.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "StMixedSoftPion.h"
/* **************************************************
 *
 * Event class used for mixed event buffer, stripped down
 * to minimum information neede to reconstruct the helix
 * and basic track information. Currently include:
 * 1) primVtx
 * 2) B-Field
 *
 * **************************************************
 *
 *  Initial Authors:
 *         ** Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

class StMixedD0;
class StMixedD0Pion;
//class StMixedSoftPion;
class StDstarMixedEvent
{
public:
   StDstarMixedEvent();
   StDstarMixedEvent(StDstarMixedEvent*);
   StDstarMixedEvent(StThreeVectorF const& vertexPos, float B, float weight = 1);
   ~StDstarMixedEvent()
   {
      ;
   };
   void addSoftPion(StMixedSoftPion const &);
   void addD0(StMixedD0 const &);
   void setPos(float,float,float);
   void setField(float);
   int getNoD0s() const;
   int getNoPions() const;
   StMixedSoftPion const&  pionAt(int) const;
   StMixedD0 const&  D0At(int) const;
   StThreeVectorF const & vertex() const;
   double const field() const;
   float const weight() const;
private:
   StThreeVectorF mVtx;
   float mBField;
   float mWeight;
   std::vector <StMixedD0 > mD0s;
   std::vector <StMixedSoftPion > mPions;
};
inline void StDstarMixedEvent::setPos(float const vx, float const vy, float const vz)
{
   mVtx = StThreeVectorF(vx, vy, vz);
}
inline void StDstarMixedEvent::setField(float const field)
{
   mBField = field;
}
inline int StDstarMixedEvent::getNoPions() const
{
   return mPions.size();
}
inline int StDstarMixedEvent::getNoD0s() const
{
   return mD0s.size();
}
inline StMixedSoftPion const& StDstarMixedEvent::pionAt(int counter) const
{
   return mPions[counter];
}
inline StMixedD0 const& StDstarMixedEvent::D0At(int counter) const
{
   return mD0s[counter];
}
inline StThreeVectorF const & StDstarMixedEvent::vertex() const
{
   return mVtx;
}
inline double const StDstarMixedEvent::field() const
{
   return mBField;
}
inline float const StDstarMixedEvent::weight() const
{
   return mWeight;
}

#endif
