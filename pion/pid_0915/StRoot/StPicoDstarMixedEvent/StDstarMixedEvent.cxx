//not include the eventplane

#include "StDstarMixedEvent.h"
#include "StMixedD0.h"
#include "StMixedSoftPion.h"
//#include "StEventPlane/StEventPlane.h"
#include "TVector2.h"
#include <limits>
StDstarMixedEvent::StDstarMixedEvent() :  mVtx(StThreeVectorF()),
   mBField(std::numeric_limits<float>::quiet_NaN())
{
}
StDstarMixedEvent::StDstarMixedEvent(StDstarMixedEvent *t) : mVtx(t->mVtx), mBField(t->mBField),mWeight(t->mWeight)
{
}
StDstarMixedEvent::StDstarMixedEvent(StThreeVectorF const& vtx, float b, float weight) :
  mVtx(vtx), mBField(b), mWeight(weight)
{
}
void StDstarMixedEvent::addD0(StMixedD0 const& D0)
{
   mD0s.push_back(D0);
   return;
}
void StDstarMixedEvent::addSoftPion(StMixedSoftPion const& pion)
{
   mPions.push_back(pion);
   return;
}