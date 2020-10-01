// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef DETECTOR_DCS_DCSPROCESSOR_H_
#define DETECTOR_DCS_DCSPROCESSOR_H_

#include <memory>
#include <Rtypes.h>
#include <unordered_map>
#include <deque>
#include "Framework/Logger.h"
#include "DetectorsDCS/DataPointCompositeObject.h"
#include "DetectorsDCS/DataPointIdentifier.h"
#include "DetectorsDCS/DataPointValue.h"
#include "DetectorsDCS/DeliveryType.h"

/// @brief Class to process DCS data points

namespace o2
{
namespace dcs
{

class DCSProcessor
{

  using Ints = std::vector<int>;
  using Chars = std::vector<char>;
  using Doubles = std::vector<double>;
  using Binaries = std::array<uint64_t, 7>;
  using Strings = std::array<char, 56>;

  using DQChars = std::deque<char>;
  using DQInts = std::deque<int>;
  using DQDoubles = std::deque<double>;
  using DQUInts = std::deque<uint32_t>;
  using DQBools = std::deque<bool>;
  using DQStrings = std::deque<Strings>;
  using DQTimes = std::deque<uint32_t>;
  using DQBinaries = std::deque<Binaries>;

  using DPID = o2::dcs::DataPointIdentifier;
  using DPVAL = o2::dcs::DataPointValue;
  using DPCOM = o2::dcs::DataPointCompositeObject;

 public:
  DCSProcessor() = default;
  ~DCSProcessor() = default;

  void init(const std::vector<DPID>& aliaseschars, const std::vector<DPID>& aliasesints,
            const std::vector<DPID>& aliasesdoubles, const std::vector<DPID>& aliasesUints,
            const std::vector<DPID>& aliasesbools, const std::vector<DPID>& aliasesstrings,
            const std::vector<DPID>& aliasestimes, const std::vector<DPID>& aliasesbinaries);

  void init(const std::vector<DPID>& aliases);

  int process(const std::unordered_map<DPID, DPVAL>& map);

  void processAlias(DPID& alias, DeliveryType type, const std::unordered_map<DPID, DPVAL>& map, std::unordered_map<DPID, DPVAL>::const_iterator& it);

  void processChars();
  void processInts();
  void processDoubles();
  void processUInts();
  void processBools();
  void processStrings();
  void processTimes();
  void processBinaries();
  uint64_t processFlag(uint64_t flag, const char* alias);

  void doSimpleMovingAverage(int nelements, std::deque<int>& vect, float& avg, bool& isSMA);

  DQChars& getVectorForAliasChar(const DPID& id) { return mDpscharsmap[id]; }
  DQInts& getVectorForAliasInt(const DPID& id) { return mDpsintsmap[id]; }
  DQDoubles& getVectorForAliasDouble(const DPID& id) { return mDpsdoublesmap[id]; }
  DQUInts& getVectorForAliasUInt(const DPID& id) { return mDpsUintsmap[id]; }
  DQBools& getVectorForAliasBool(const DPID& id) { return mDpsboolsmap[id]; }
  DQStrings& getVectorForAliasString(const DPID& id) { return mDpsstringsmap[id]; }
  DQTimes& getVectorForAliasTime(const DPID& id) { return mDpstimesmap[id]; }
  DQBinaries& getVectorForAliasBinary(const DPID& id) { return mDpsbinariesmap[id]; }

 private:
  float mAvgTestInt0 = 0; // moving average for DP named TestInt0
  std::unordered_map<DPID, DQChars> mDpscharsmap;
  std::unordered_map<DPID, DQInts> mDpsintsmap;
  std::unordered_map<DPID, DQDoubles> mDpsdoublesmap;
  std::unordered_map<DPID, DQUInts> mDpsUintsmap;
  std::unordered_map<DPID, DQBools> mDpsboolsmap;
  std::unordered_map<DPID, DQStrings> mDpsstringsmap;
  std::unordered_map<DPID, DQTimes> mDpstimesmap;
  std::unordered_map<DPID, DQBinaries> mDpsbinariesmap;
  std::vector<DPID> mAliaseschars;
  std::vector<DPID> mAliasesints;
  std::vector<DPID> mAliasesdoubles;
  std::vector<DPID> mAliasesUints;
  std::vector<DPID> mAliasesbools;
  std::vector<DPID> mAliasesstrings;
  std::vector<DPID> mAliasestimes;
  std::vector<DPID> mAliasesbinaries;
  std::vector<uint64_t> mLatestTimestampchars;
  std::vector<uint64_t> mLatestTimestampints;
  std::vector<uint64_t> mLatestTimestampdoubles;
  std::vector<uint64_t> mLatestTimestampUints;
  std::vector<uint64_t> mLatestTimestampbools;
  std::vector<uint64_t> mLatestTimestampstrings;
  std::vector<uint64_t> mLatestTimestamptimes;
  std::vector<uint64_t> mLatestTimestampbinaries;

  ClassDefNV(DCSProcessor, 0);
};

} // namespace dcs
} // namespace o2

#endif
