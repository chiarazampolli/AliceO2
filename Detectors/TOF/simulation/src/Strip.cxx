// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//
//  Strip.cxx: structure to store the TOF hits and digits in strips - useful
// for clusterization purposes
//  ALICEO2
//
#include <cstring>
#include <tuple>

#include <TMath.h>
#include <TObjArray.h>

#include "TOFSimulation/Strip.h"

using namespace o2::tof;

ClassImp(Strip);

//_______________________________________________________________________
Strip::Strip(Int_t index)
  : mStripIndex(index)
{
}
//_______________________________________________________________________
Strip::Strip(const Strip& ref) = default;

//_______________________________________________________________________
void Strip::InsertHit(const HitType* h)
{
  int stripIndex = GetStripIndex(h);
  if (stripIndex != mStripIndex) {
    LOG(ERROR) << "Strip Index of of requested hit (" << stripIndex << ") does not match with current strip (" << nStripIndex << ")" << endl;
    return;
  }
  mHits.push_back(h);
}

//_______________________________________________________________________
const HitType* Strip::GetHitAt(Int_t i) const
{
  if (i < mHits.size()) {
    return mHits[i];
  }
  return nullptr;
}

//_______________________________________________________________________
void Strip::Clear() { ClearHits(); }
//_______________________________________________________________________
Int_t Strip::GetStripIndex(const HitType* hit)
{
  // finds the strip index given the hit
  Float_t pos[3] = { hit.GetX(), hit.GetY(), hit.GetZ() };
  Int_t detInd[5];
  Geo::getDetID(pos, detInd);
  Int_t channelID = Geo::getIndex(detInd);
  return channelID/Geo::NPADS;

}
//_______________________________________________________________________
Int_t Strip::addDigit(Double_t time, Int_t channel, Int_t tdc, Int_t tot, Int_t bc, Int_t lbl)
{

  // return the MC label. We pass it also as argument, but it can change in
  // case the digit was merged
  
  auto key = Digit::getOrderingKey(channel, bc, tdc);
  auto dig = findDigit(key);
  if (dig) {
    lbl = dig->getLabel();
    dig->addCharge(charge, lbl);
  } else {
    auto digIter = mDigits.emplace(
				   //std::make_pair(key, Digit(static_cast<UShort_t>(mChipIndex), roframe, row, col, charge, timestamp)));
				   std::make_pair(key, Digit(time, channel, tdc, tot, bc, lbl)));
    auto pair = digIter.first;
    dig = &(pair->second);
    dig->setLabel(0, lbl);
  }
  return lbl;
}

//______________________________________________________________________



