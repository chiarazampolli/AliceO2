// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TOF_DIGITIZER_H_
#define ALICEO2_TOF_DIGITIZER_H_

#include "TOFBase/Geo.h"
#include "TOFBase/Digit.h"
#include "TOFSimulation/Detector.h"
#include "TOFSimulation/Strip.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TOFSimulation/MCLabel.h"

namespace o2
{
namespace tof
{
class Digitizer
{
 public:
  Digitizer(Int_t mode = 0) : mMode(mode), mReadoutWindowCurrent(0) { init(); };
  ~Digitizer() = default;

  void init();

  void process(const std::vector<HitType>* hits, std::vector<Digit>* digits);

  Float_t getShowerTimeSmeared(Float_t time, Float_t charge);
  Float_t getDigitTimeSmeared(Float_t time, Float_t x, Float_t z, Float_t charge);
  Float_t getCharge(Float_t eDep);
  Bool_t isFired(Float_t x, Float_t z, Float_t charge);
  Float_t getEffX(Float_t x);
  Float_t getEffZ(Float_t z);
  Float_t getFractionOfCharge(Float_t x, Float_t z);

  Int_t getCurrentReadoutWindow() const { return mReadoutWindowCurrent; }
  void setCurrentReadoutWindow(Double_t value) { mReadoutWindowCurrent = value; }
  Float_t getTimeLastHit(Int_t idigit) const { return 0; }
  Float_t getTotLastHit(Int_t idigit) const { return 0; }
  Int_t getXshift(Int_t idigit) const { return 0; }
  Int_t getZshift(Int_t idigit) const { return 0; }
  void setEventTime(double value) { mEventTime = value; }
  void setEventID(Int_t id) { mEventID = id; }
  void setSrcID(Int_t id) { mSrcID = id; }

  void setMCTruthContainer(o2::dataformats::MCTruthContainer<o2::MCCompLabel>* truthcontainer){
    mMCTrueContainer = truthcontainer;
  }

  void initParameters();
  void printParameters();

  void test(const char* geo = "O2geometry.root");
  void testFromHits(const char* geo = "O2geometry.root", const char* hits = "AliceO2_TGeant3.tof.mc_10_event.root");
  void   fillOutputContainer(std::vector<Digit>* digits);

  void setContinuous(bool val) {mContinuous=val;}

 private:
  // parameters
  Int_t mMode;
  Float_t mBound1;
  Float_t mBound2;
  Float_t mBound3;
  Float_t mBound4;
  Float_t mTOFresolution;
  Float_t mShowerResolution;
  Float_t mDigitResolution;
  Float_t mTimeSlope;
  Float_t mTimeDelay;
  Float_t mTimeDelayCorr;
  Float_t mTimeWalkeSlope;
  Float_t mEffCenter;
  Float_t mEffBoundary1;
  Float_t mEffBoundary2;
  Float_t mEffBoundary3;

  // info TOF timewindow
  Int_t mReadoutWindowCurrent;
  Double_t mEventTime;
  Int_t mEventID = 0;
  Int_t mSrcID = 0;

  bool mContinuous=false;

  // digit info
  std::vector<Digit>* mDigits;
  o2::dataformats::MCTruthContainer<o2::tof::MCLabel> mMCTruthContainerA;
  o2::dataformats::MCTruthContainer<o2::tof::MCLabel> mMCTruthContainerB;
  o2::dataformats::MCTruthContainer<o2::tof::MCLabel>* mMCTruthContainerCurrent = &mMCTruthContainerA; ///< Array for MCTruth information associated to digits in mDigitsArrray. 
  o2::dataformats::MCTruthContainer<o2::tof::MCLabel>* mMCTruthContainerNext = &mMCTruthContainerB; ///< Array for MCTruth information associated to digits in mDigitsArrray. 
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mMCTrueContainer;

  // array of strips to store the digits per strip (one for the current readout window, one for the next one)
  std::vector<Strip> mStripsA;
  std::vector<Strip> mStripsB;
  std::vector<Strip>* mStripsCurrent = &(mStripsA);
  std::vector<Strip>* mStripsNext = &(mStripsB);

  
  Int_t processHit(const HitType& hit, Double_t event_time);
  void addDigit(Int_t channel, UInt_t istrip, Float_t time, Float_t x, Float_t z, Float_t charge, Int_t iX, Int_t iZ, Int_t padZfired,
                Int_t trackID);

  bool isMergable(Digit digit1, Digit digit2)
  {
    if (digit1.getChannel() != digit2.getChannel()) {
      return false;
    }

    if (digit1.getBC() != digit2.getBC()) {
      return false;
    }

    // Check if the difference is larger than the TDC dead time
    if (std::abs(digit1.getTDC() - digit2.getTDC()) > Geo::DEADTIMETDC) {
      return false;
    }
    return true;
  }

  ClassDefNV(Digitizer, 1);
};
}
}
#endif
