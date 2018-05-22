// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HitReader.h
/// \brief Definition of the TOF hit reader

#ifndef ALICEO2_TOF_HITREADER_H
#define ALICEO2_TOF_HITREADER_H

#include "TOFBase/Digit.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2
{
namespace TOF
{
/// \class HitReader
/// \brief HitReader class for TOF
///
class HitReader
{
  using Label = o2::MCCompLabel;

  HitReader() = default;
  HitReader(const HitReader& cluster) = delete;
  virtual ~HitReader() = default;

  HitReader& operator=(const HitReader& src) = delete;

  virtual void init() = 0;
  //
 protected:
  //
};

/// \class DigitHitReader
/// \brief DigitHitReader class for TOF. Feeds the MC digits to the Cluster Finder
///
class DigitHitReader : public HitReader
{
 public:
  DigitHitReader() = default;
  void setDigitArray(const std::vector<o2::TOF::Digit>* a)
  {
    mDigitArray = a;
    mIdx = 0;
  }

  void init() override
  {
    mIdx = 0;
    mLastDigit = nullptr;
  }

 private:
  const std::vector<o2::ITSMFT::Digit>* mDigitArray = nullptr;
  const Digit* mLastDigit = nullptr;
  Int_t mIdx = 0;
};

/// \class RawHitReader
/// \brief RawHitReader class for TOF. Feeds raw data to the Cluster Finder
///
class RawHitReader : public HitReader
{
 public:
  Bool_t getNextChipData(ChipPixelData& chipData) override;
};

} // namespace ITSMFT
} // namespace o2

#endif /* ALICEO2_ITS_PIXELREADER_H */
