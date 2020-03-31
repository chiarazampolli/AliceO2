// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef TOF_CHANNEL_CALIBRATOR_H_
#define TOF_CHANNEL_CALIBRATOR_H_

#include "DetectorsCalibration/TimeSlotCalibration.h"
#include "DetectorsCalibration/TimeSlot.h"
#include "DataFormatsTOF/CalibInfoTOF.h"
#include "TOFCalibration/CalibTOFapi.h"
#include "DataFormatsTOF/CalibLHCphaseTOF.h"
#include "TOFBase/Geo.h"
#include "CCDB/CcdbObjectInfo.h"
#include <array>
#include <boost/histogram.hpp>

namespace o2
{
namespace tof
{

class TOFChannelData {

using Slot = o2::calibration::TimeSlot<o2::tof::TOFChannelData>;

public: 
  TOFChannelData() {
    LOG(INFO) << "Default c-tor, not to be used";
  }

  TOFChannelData(int nb, float r) : mNbins(nb), mRange(r)
  {
    if (r <= 0. || nb < 1) {
      throw std::runtime_error("Wrong initialization of the histogram");
    }
    mV2Bin = mNbins / (2 * mRange);
    mHisto = boost::histogram::make_histogram(boost::histogram::axis::regular<>(mNbins, -mRange, mRange, "t-text"),
					      boost::histogram::axis::regular<>(o2::tof::Geo::NCHANNELS, -0.5, o2::tof::Geo::NCHANNELS-0.5, "channel index")); // bin along channel axis is centered in the channel index
  }

  ~TOFChannelData() = default;

  void print() const;
  void fill(const gsl::span<const o2::dataformats::CalibInfoTOF> data);
  void merge(const TOFChannelData* prev);
  int findBin(float v) const;
  float integral(int ch, int binmin, int binmax) const;
  bool hasEnoughData(const Slot& slot, int minEntries) const;

  float getRange() const { return mRange; }
  void setRange(float r) { mRange = r; }

  float getNbins() const { return mNbins; }
  void setNbins(float nb) { mNbins = nb; }

  boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double, boost::use_default, boost::use_default, boost::use_default> >, boost::histogram::unlimited_storage<std::allocator<char> > > getHisto() const { return mHisto; }
    
 private:
  float mRange = 24400;
  int mNbins = 1000;
  float mV2Bin;
  // do I really have to initialize it like below?
  boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double, boost::use_default, boost::use_default, boost::use_default> >, boost::histogram::unlimited_storage<std::allocator<char> > > mHisto;

  ClassDefNV(TOFChannelData, 1);
};

class TOFChannelCalibrator : public o2::calibration::TimeSlotCalibration<o2::dataformats::CalibInfoTOF, o2::tof::TOFChannelData>
{
  using TFType = uint64_t;
  using Slot = o2::calibration::TimeSlot<o2::tof::TOFChannelData>;
  using CalibTOFapi = o2::tof::CalibTOFapi;
  using TimeSlewing = o2::dataformats::CalibTimeSlewingParamTOF;
  using CcdbObjectInfo = o2::ccdb::CcdbObjectInfo;
  using CcdbObjectInfoVector = std::vector<CcdbObjectInfo>;
  using TimeSlewingVector = std::vector<TimeSlewing>;

 public:
  TOFChannelCalibrator(int minEnt = 500, int nb = 1000, float r = 24400, const std::string path = "http://ccdb-test.cern.ch:8080") : mMinEntries(minEnt), mNBins(nb), mRange(r) { mCalibTOFapi.setURL(path); }

  bool hasEnoughData(const Slot& slot) const final;
  void initOutput() final;
  void finalizeSlot(Slot& slot) final;
  Slot& emplaceNewSlot(bool front, TFType tstart, TFType tend) final;

  const TimeSlewingVector& getTimeSlewingVector() const { return mTimeSlewingVector; }
  const CcdbObjectInfoVector& getTimeSlewingInfoVector() const { return mInfoVector; }
  CcdbObjectInfoVector& getTimeSlewingInfoVector() { return mInfoVector; }

 private:
  int mMinEntries = 0;  // min number of entries to calibrate the TimeSlot
  int mNBins = 0;  // bins of the histogram with the t-text per channel
  float mRange = 0.;  // range of the histogram with the t-text per channel
  CalibTOFapi mCalibTOFapi; 
  CcdbObjectInfoVector mInfoVector; // vector of CCDB Infos , each element is filled with the CCDB description of the accompanying TimeSlewing object
  TimeSlewingVector mTimeSlewingVector;   // vector of TimeSlewing, each element is filled in "process"
                                          // when we finalize one slot (multiple can be finalized
                                          // during the same "process", which is why we have a vector).
                                          // Each element is to be considered the output of the device,
                                          // and will go to the CCDB. Note that for the channel offset
                                          // we still fill the TimeSlewing object

  ClassDefOverride(TOFChannelCalibrator, 1);
};

} // end namespace tof
} // end namespace o2

#endif /* TOF_CHANNEL_CALIBRATOR_H_ */
