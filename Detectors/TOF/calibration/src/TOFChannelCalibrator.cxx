// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TOFCalibration/TOFChannelCalibrator.h"
#include "Framework/Logger.h"
#include "MathUtils/MathBase.h"
#include "CommonUtils/MemFileHelper.h"
#include "CCDB/CcdbApi.h"
#include "DetectorsCalibration/Utils.h"
#include <boost/histogram/ostream.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <iostream>
#include <sstream>

namespace o2
{
namespace tof
{

using Slot = o2::calibration::TimeSlot<o2::tof::TOFChannelData>;
using TimeSlewing = o2::dataformats::CalibTimeSlewingParamTOF;
using clbUtils = o2::calibration::Utils;
using o2::math_utils::math_base::fitGaus;

//_____________________________________________
TOFChannelData::TOFChannelData()
{
  LOG(INFO) << "Default c-tor, not to be used";
}

//_____________________________________________
void TOFChannelData::fill(const gsl::span<const o2::dataformats::CalibInfoTOF> data)
{
  // fill container
  for (int i = data.size(); i--;) {
    auto dt = data[i].getDeltaTimePi();
    auto ch = data[i].getTOFChIndex();
    histo.fill(dt, ch);
  }
}

//_____________________________________________
void TOFChannelData::merge(const LHCClockDataHisto* prev)
{
  // merge data of 2 slots
  histo += prev->histo;
}
  
//_____________________________________________
bool TOFChannelData::hasEnoughData(const Slot& slot) const
{
  // true if all channels can be fitted --> have enough statistics
  for (ich = 0; ich < o2::tof::Geo::NCHANNELS; ich++) {
    // make the slice of the 2D histogram so that we have the 1D of the current channel
    auto hch = algorithm::reduce(c->histo, algorithm::shrink(1, float(ich), float(ich)));
    if (hch.sum() < mMinEntries) return false;
  }
  return true;
}
  
//_____________________________________________
void TOFChannelData::print() const
{
  LOG(INFO) << "Printing histogram:";
  std::ostringstream os;
  os << histo;
  LOG(INFO) << "Number of entries in histogram: " << histo.sum();
}

//_____________________________________________
int TOFChannelData::findBin(float v) const
{
  // find the bin along the x-axis (with t-texp) where the value "v" is; this does not depend on the channel (axis 1)
  for (auto&& x : indexed(histo)) {
    const auto b0 = x.bin(0); // we take the bin limits along axis 0
    if (b0.lower() <= v && b0.upper() > v) return x.index(0);
  }
  return -1; // if we cannot find the bin  
}

//_____________________________________________
float TOFChannelData::integral(int ch, int binmin, int binmax) const
{
  // calculates the integral along one channel only, in [binmin, binmax]
  auto hch = algorithm::reduce(histo, algorithm::shrink(1, float(ich), float(ich)), algorithm::slice(0, binmin, binmax));
  return hch.sum();
}
  
//===================================================================

//_____________________________________________
void TOFChannelCalibrator::initOutput()
{
  // Here we initialize the vector of our output objects
  mInfoVector.clear();
  mTimeSlewingVector.clear();
  return;
}

//_____________________________________________
void TOFChannelCalibrator::finalizeSlot(Slot& slot)
{
  // Extract results for the single slot
  o2::tof::TOFChannelData* c = slot.getContainer();
  LOG(INFO) << "Finalize slot " << slot.getTFStart() << " <= TF <= " << slot.getTFEnd() << " with "
	    << c->getEntries() << " entries";

  // for the CCDB entry
  std::map<std::string, std::string> md;
  TimeSlewing ts;
  
  for (ich = 0; ich < o2::tof::Geo::NCHANNELS; ich++) {
    // make the slice of the 2D histogram so that we have the 1D of the current channel
    auto hch = algorithm::reduce(c->histo, algorithm::shrink(1, float(ich), float(ich)));
    std::vector<float> fitValues;
    std::vector<float> histoValues;
    for (auto x : indexed(hch)) {
      histoValues.push_back(x.get());
    }
    int fitres = fitGaus(c->nbins, histoValues.data(), -(c->range), c->range, fitValues);
    if (fitres >= 0) {
      LOG(INFO) << "Channel " << ich << " :: Fit result " << fitres << " Mean = " << fitValues[1] << " Sigma = " << fitValues[2];
    } else {
      LOG(ERROR) << "Channel " << ich << " :: Fit failed with result = " << fitres;
    }
    float fractionUnderPeak;
    float intmin = fitValues[1] - 5 * fitValues[2]; // mean - 5*sigma
    float intmax = fitValues[1] + 5 * fitValues[2]; // mean + 5*sigma
    float intmin2;
    float intmax2;
    
    // if peak is at the border of our bunch-crossing window (-12.5:12.5 ns)
    // continue to extrapolate gaussian integral from the other border
    float addduetoperiodicity = 0;
    if (intmin < -12500) { // at left border
      intmin2 = intmin + 25000;
      intmin = -12500;
      intmax2 = 12500;
      if (intmin2 > intmax) {
	int binmin2 = c->findBin(intmin2);
	int binmax2 = c->findBin(intmax2);
	addduetoperiodicity = c->integral(ich, binmin2, binmax2);
      }
    } else if (intmax > 12500) { // at right border
      intmax2 = intmax - 25000;
      intmax = 12500;
      intmin2 = -12500;
      if (intmax2 < intmin) {
	int binmin2 = c->findBin(intmin2);
	int binmax2 = c->findBin(intmax2);
	addduetoperiodicity = c->integral(ich, binmin2, binmax2);
      }
    }
    
    int binmin = c->findBin(intmin);
    int binmax = c->findBin(intmax);

    if (binmin < 1)
      binmin = 1; // avoid to take the underflow bin (can happen in case the sigma is too large)
    if (binmax > c->nbins)
      binmax = c->nbins; // avoid to take the overflow bin (can happen in case the sigma is too large)

    float fractionUnderPeak = (c->integral(ch, binmin, binmax) + addduetoperiodicity) / c->integral(ch, 1, c->nbins());
    // now we need to store the results in the TimeSlewingObject
    ts.setFractionUnderPeak(ich / o2::tof::Geo::NPADSXSECTOR, ich % o2::tof::Geo::NPADSXSECTOR, fractionUnderPeak);
    ts.setSigmaPeak(ich / o2::tof::Geo::NPADSXSECTOR, ich % o2::tof::Geo::NPADSXSECTOR, abs(fitValues[2]));
    ts.addTimeSlewingInfo(ich, 0, fitValues[1]);
  }
  auto clName = o2::utils::MemFileHelper::getClassName(ts);
  auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
  mInfoVector.emplace_back("TOF/ChannelCalib", clName, flName, md, slot.getTFStart(), 99999999999999);
  mTimeSlewingVector.emplace_back(l);

  slot.print();
}

//_____________________________________________
Slot& TOFChannelCalibrator::emplaceNewSlot(bool front, TFType tstart, TFType tend)
{
  auto& cont = getSlots();
  auto& slot = front ? cont.emplace_front(tstart, tend) : cont.emplace_back(tstart, tend);
  slot.setContainer(std::make_unique<TOFChannelData>(mNBins, mRange));
  return slot;
}

} // end namespace tof
} // end namespace o2
  
