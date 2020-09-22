// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <DetectorsDCS/DCSProcessor.h>
#include "Rtypes.h"
#include <deque>
#include <string.h>

using namespace o2::dcs;

using DeliveryType = o2::dcs::DeliveryType;
using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;

//ClassImp(o2::dcs::DCSProcessor);

void DCSProcessor::init(std::vector<DPID> aliaseschars, std::vector<DPID> aliasesints, std::vector<DPID> aliasesdoubles,
                        std::vector<DPID> aliasesUints, std::vector<DPID> aliasesbools, std::vector<DPID> aliasesstrings,
                        std::vector<DPID> aliasestimes, std::vector<DPID> aliasesbinaries)
{

  // init from separate vectors of aliases (one per data point type)

  // chars
  for (auto it = std::begin(aliaseschars); it != std::end(aliaseschars); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_CHAR) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a char";
    }
    mAliaseschars.emplace_back((*it).get_alias(), DeliveryType::RAW_CHAR);
  }

  // ints
  for (auto it = std::begin(aliasesints); it != std::end(aliasesints); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_INT) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a int";
    }
    mAliasesints.emplace_back((*it).get_alias(), DeliveryType::RAW_INT);
  }

  // doubles
  for (auto it = std::begin(aliasesdoubles); it != std::end(aliasesdoubles); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_DOUBLE) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a double";
    }
    mAliasesdoubles.emplace_back((*it).get_alias(), DeliveryType::RAW_DOUBLE);
  }

  // uints
  for (auto it = std::begin(aliasesUints); it != std::end(aliasesUints); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_UINT) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a uint";
    }
    mAliasesUints.emplace_back((*it).get_alias(), DeliveryType::RAW_UINT);
  }

  // bools
  for (auto it = std::begin(aliasesbools); it != std::end(aliasesbools); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_BOOL) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a bool";
    }
    mAliasesbools.emplace_back((*it).get_alias(), DeliveryType::RAW_BOOL);
  }

  // strings
  for (auto it = std::begin(aliasesstrings); it != std::end(aliasesstrings); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_STRING) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a string";
    }
    mAliasesstrings.emplace_back((*it).get_alias(), DeliveryType::RAW_STRING);
  }

  // times
  for (auto it = std::begin(aliasestimes); it != std::end(aliasestimes); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_TIME) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a time";
    }
    mAliasestimes.emplace_back((*it).get_alias(), DeliveryType::RAW_TIME);
  }

  // binaries
  for (auto it = std::begin(aliasesbinaries); it != std::end(aliasesbinaries); ++it) {
    if ((*it).get_type() != DeliveryType::RAW_BINARY) {
      LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a binary";
    }
    mAliasesbinaries.emplace_back((*it).get_alias(), DeliveryType::RAW_BINARY);
  }

  mLatestTimestampchars.resize(aliaseschars.size());
  mLatestTimestampints.resize(aliasesints.size());
  mLatestTimestampdoubles.resize(aliasesdoubles.size());
  mLatestTimestampUints.resize(aliasesUints.size());
  mLatestTimestampbools.resize(aliasesbools.size());
  mLatestTimestampstrings.resize(aliasesstrings.size());
  mLatestTimestamptimes.resize(aliasestimes.size());
  mLatestTimestampbinaries.resize(aliasesbinaries.size());

  for (auto i = 0; i < aliaseschars.size(); i++) {
    mLatestTimestampchars[i] = 0;
  }
  for (auto i = 0; i < aliasesints.size(); i++) {
    mLatestTimestampints[i] = 0;
  }
  for (auto i = 0; i < aliasesdoubles.size(); i++) {
    mLatestTimestampdoubles[i] = 0;
  }
  for (auto i = 0; i < aliasesUints.size(); i++) {
    mLatestTimestampUints[i] = 0;
  }
  for (auto i = 0; i < aliasesbools.size(); i++) {
    mLatestTimestampbools[i] = 0;
  }
  for (auto i = 0; i < aliasesstrings.size(); i++) {
    mLatestTimestampstrings[i] = 0;
  }
  for (auto i = 0; i < aliasestimes.size(); i++) {
    mLatestTimestamptimes[i] = 0;
  }
  for (auto i = 0; i < aliasesbinaries.size(); i++) {
    mLatestTimestampbinaries[i] = 0;
  }
}

//______________________________________________________________________

void DCSProcessor::init(std::vector<DPID> aliases)
{

  int nchars = 0, nints = 0, ndoubles = 0, nUints = 0,
      nbools = 0, nstrings = 0, ntimes = 0, nbinaries = 0;
  for (auto it = std::begin(aliases); it != std::end(aliases); ++it) {
    if ((*it).get_type() == DeliveryType::RAW_CHAR) {
      mAliaseschars.emplace_back((*it).get_alias(), DeliveryType::RAW_CHAR);
      nchars++;
    }
    if ((*it).get_type() == DeliveryType::RAW_INT) {
      mAliasesints.emplace_back((*it).get_alias(), DeliveryType::RAW_INT);
      nints++;
    }
    if ((*it).get_type() == DeliveryType::RAW_DOUBLE) {
      mAliasesdoubles.emplace_back((*it).get_alias(), DeliveryType::RAW_DOUBLE);
      ndoubles++;
    }
    if ((*it).get_type() == DeliveryType::RAW_UINT) {
      mAliasesUints.emplace_back((*it).get_alias(), DeliveryType::RAW_UINT);
      nUints++;
    }
    if ((*it).get_type() == DeliveryType::RAW_BOOL) {
      mAliasesbools.emplace_back((*it).get_alias(), DeliveryType::RAW_BOOL);
      nbools++;
    }
    if ((*it).get_type() == DeliveryType::RAW_STRING) {
      mAliasesstrings.emplace_back((*it).get_alias(), DeliveryType::RAW_STRING);
      nstrings++;
    }
    if ((*it).get_type() == DeliveryType::RAW_TIME) {
      mAliasestimes.emplace_back((*it).get_alias(), DeliveryType::RAW_TIME);
      ntimes++;
    }
    if ((*it).get_type() == DeliveryType::RAW_BINARY) {
      mAliasesbinaries.emplace_back((*it).get_alias(), DeliveryType::RAW_BINARY);
      nbinaries++;
    }
  }

  mLatestTimestampchars.resize(nchars);
  mLatestTimestampints.resize(nints);
  mLatestTimestampdoubles.resize(ndoubles);
  mLatestTimestampUints.resize(nUints);
  mLatestTimestampbools.resize(nbools);
  mLatestTimestampstrings.resize(nstrings);
  mLatestTimestamptimes.resize(ntimes);
  mLatestTimestampbinaries.resize(nbinaries);

  for (auto i = 0; i < nchars; i++) {
    mLatestTimestampchars[i] = 0;
  }
  for (auto i = 0; i < nints; i++) {
    mLatestTimestampints[i] = 0;
  }
  for (auto i = 0; i < ndoubles; i++) {
    mLatestTimestampdoubles[i] = 0;
  }
  for (auto i = 0; i < nUints; i++) {
    mLatestTimestampUints[i] = 0;
  }
  for (auto i = 0; i < nbools; i++) {
    mLatestTimestampbools[i] = 0;
  }
  for (auto i = 0; i < nstrings; i++) {
    mLatestTimestampstrings[i] = 0;
  }
  for (auto i = 0; i < ntimes; i++) {
    mLatestTimestamptimes[i] = 0;
  }
  for (auto i = 0; i < nbinaries; i++) {
    mLatestTimestampbinaries[i] = 0;
  }
}

//__________________________________________________________________

int DCSProcessor::process(const std::unordered_map<DPID, DPVAL>& map)
{

  // process function to do "something" with the DCS map that is passed

  // first, we need to check if there are the Data Points that we need

  auto foundChars = 0, foundInts = 0, foundDoubles = 0, foundUInts = 0,
       foundBools = 0, foundStrings = 0, foundTimes = 0, foundBinaries = 0;

  // char type
  auto s = mAliaseschars.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliaseschars.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliaseschars[i], DeliveryType::RAW_CHAR, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundChars++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliaseschars[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpscharsmap[mAliaseschars[i]].size() = " << mDpscharsmap[mAliaseschars[i]].size();
        if (mDpscharsmap[mAliaseschars[i]].size() == 0 || (i > 0 && etime != mLatestTimestampchars[i])) {
          mDpscharsmap[mAliaseschars[i]].push_back(val.payload_pt1);
          mLatestTimestampchars[i] = etime;
        }
      }
    }
    processChars();
  }

  // int type
  s = mAliasesints.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasesints.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasesints[i], DeliveryType::RAW_INT, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundInts++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasesints[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpsintsmap[mAliasesints[i]].size() = " << mDpsintsmap[mAliasesints[i]].size();
        if (mDpsintsmap[mAliasesints[i]].size() == 0 || (etime != mLatestTimestampints[i])) {
          mDpsintsmap[mAliasesints[i]].push_back(val.payload_pt1);
          mLatestTimestampints[i] = etime;
        }
      }
    }
    processInts();
  }

  // double type
  s = mAliasesdoubles.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasesdoubles.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasesdoubles[i], DeliveryType::RAW_DOUBLE, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundDoubles++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasesdoubles[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpsdoublemap[mAliasesdoubles[i]].size() = " << mDpsdoublesmap[mAliasesdoubles[i]].size();
        if (mDpsdoublesmap[mAliasesdoubles[i]].size() == 0 || (etime != mLatestTimestampdoubles[i])) {
          mDpsdoublesmap[mAliasesdoubles[i]].push_back(val.payload_pt1);
          mLatestTimestampdoubles[i] = etime;
        }
      }
    }
    processDoubles();
  }

  // UInt type
  s = mAliasesUints.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasesUints.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasesUints[i], DeliveryType::RAW_UINT, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundUInts++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasesUints[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpsUintsmap[mAliasesUints[i]].size() = " << mDpsUintsmap[mAliasesUints[i]].size();
        if (mDpsUintsmap[mAliasesUints[i]].size() == 0 || (etime != mLatestTimestampUints[i])) {
          mDpsUintsmap[mAliasesUints[i]].push_back(val.payload_pt1);
          mLatestTimestampUints[i] = etime;
        }
      }
    }
    processUInts();
  }

  // Bool type
  s = mAliasesbools.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasesbools.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasesbools[i], DeliveryType::RAW_BOOL, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundBools++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasesbools[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpsboolmap[mAliasesbools[i]].size() = " << mDpsboolsmap[mAliasesbools[i]].size();
        if (mDpsboolsmap[mAliasesbools[i]].size() == 0 || (etime != mLatestTimestampbools[i])) {
          mDpsboolsmap[mAliasesbools[i]].push_back(val.payload_pt1);
          mLatestTimestampbools[i] = etime;
        }
      }
    }
    processBools();
  }

  // String type
  s = mAliasesstrings.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasesstrings.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasesstrings[i], DeliveryType::RAW_STRING, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundStrings++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasesstrings[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpsstringmap[mAliasesstrings[i]].size() = " << mDpsstringsmap[mAliasesstrings[i]].size();
        if (mDpsstringsmap[mAliasesstrings[i]].size() == 0 || (etime != mLatestTimestampstrings[i])) {
          auto& tmp = mDpsstringsmap[mAliasesstrings[i]].emplace_back();
          std::strncpy(tmp.data(), (char*)&(val.payload_pt1), 56);
          mLatestTimestampstrings[i] = etime;
        }
      }
    }
    processStrings();
  }

  // Time type
  s = mAliasestimes.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasestimes.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasestimes[i], DeliveryType::RAW_TIME, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundTimes++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasestimes[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpstimesmap[mAliasestimes[i]].size() = " << mDpstimesmap[mAliasestimes[i]].size();
        if (mDpstimesmap[mAliasestimes[i]].size() == 0 || (etime != mLatestTimestamptimes[i])) {
          mDpstimesmap[mAliasestimes[i]].push_back(val.payload_pt1);
          mLatestTimestamptimes[i] = etime;
        }
      }
    }
    processTimes();
  }

  // Binary type
  s = mAliasesbinaries.size();
  if (s > 0) {
    for (size_t i = 0; i != mAliasesbinaries.size(); ++i) {
      int count = 0;
      std::unordered_map<DPID, DPVAL>::const_iterator it;
      processAlias(mAliasesbinaries[i], DeliveryType::RAW_BINARY, map, count, it);
      if (count == 0)
        continue; // the alias was not found in the map
      foundBinaries++;
      auto& val = it->second;
      auto flags = val.get_flags();
      if (processFlag(flags, mAliasesbinaries[i].get_alias()) == 0) {
        auto etime = val.get_epoch_time();
        // fill only if new value has a timestamp different from the timestamp of the previous one
        LOG(DEBUG) << "mDpsstringmap[mAliasesbinaries[i]].size() = " << mDpsbinariesmap[mAliasesbinaries[i]].size();
        if (mDpsbinariesmap[mAliasesbinaries[i]].size() == 0 || (etime != mLatestTimestampbinaries[i])) {
          auto& tmp = mDpsbinariesmap[mAliasesbinaries[i]].emplace_back();
          memcpy(tmp.data(), &(val.payload_pt1), 7);
          mLatestTimestampbinaries[i] = etime;
        }
      }
    }
    processBinaries();
  }

  if (foundChars != mAliaseschars.size())
    LOG(INFO) << "Not all expected char-typed DPs found!";
  if (foundInts != mAliasesints.size())
    LOG(INFO) << "Not all expected int-typed DPs found!";
  if (foundDoubles != mAliasesdoubles.size())
    LOG(INFO) << "Not all expected double-typed DPs found!";
  if (foundUInts != mAliasesUints.size())
    LOG(INFO) << "Not all expected uint-typed DPs found!";
  if (foundBools != mAliasesbools.size())
    LOG(INFO) << "Not all expected bool-typed DPs found!";
  if (foundStrings != mAliasesstrings.size())
    LOG(INFO) << "Not all expected string-typed DPs found!";
  if (foundTimes != mAliasestimes.size())
    LOG(INFO) << "Not all expected time-typed DPs found!";
  if (foundBinaries != mAliasesbinaries.size())
    LOG(INFO) << "Not all expected binary-typed DPs found!";

  return 0;
}

//______________________________________________________________________

void DCSProcessor::processAlias(DPID& alias, DeliveryType type, const std::unordered_map<DPID, DPVAL>& map, int& count, std::unordered_map<DPID, DPVAL>::const_iterator& it)
{

  // processing basic checks for map: all needed aliases must be present

  LOG(INFO) << "Processing " << alias;
  it = map.find(alias);
  if (it == map.end()) {
    LOG(ERROR) << "Element " << alias << " not found " << std::endl;
    count = 0;
  }
  DeliveryType tt = alias.get_type();
  if (tt != type) {
    LOG(FATAL) << "Delivery Type of alias " << alias.get_alias() << " does not match definition in DCSProcessor (" << type << ")! Please fix";
  }
  count = 1;
  return;
}

//______________________________________________________________________

void DCSProcessor::processChars()
{

  // function to process aliases of Char type; it will just print them

  for (size_t i = 0; i != mAliaseschars.size(); ++i) {
    LOG(INFO) << "processChars: mAliaseschars[" << i << "] = " << mAliaseschars[i];
    auto& id = mAliaseschars[i];
    auto& vchar = mDpscharsmap[id];
    LOG(INFO) << "vchar size = " << vchar.size();
    for (size_t j = 0; j < vchar.size(); j++) {
      LOG(INFO) << "DP = " << mAliaseschars[i] << " , value[" << j << "] = " << vchar[j];
    }
  }
}

//______________________________________________________________________

void DCSProcessor::processInts()
{

  // function to process aliases of Int type

  for (size_t i = 0; i != mAliasesints.size(); ++i) {
    LOG(INFO) << "processInts: mAliasesints[" << i << "] = " << mAliasesints[i];
    auto& id = mAliasesints[i];
    auto& vint = getVectorForAliasInt(id);
    LOG(INFO) << "vint size = " << vint.size();
    for (size_t j = 0; j < vint.size(); j++) {
      LOG(INFO) << "DP = " << mAliasesints[i] << " , value[" << j << "] = " << vint[j];
    }
    bool isSMA = false;
    LOG(INFO) << "get alias = " << id.get_alias();
    if (strcmp(id.get_alias(), "TestInt0") == 0) {
      doSimpleMovingAverage(2, vint, mAvgTestInt0, isSMA);
      LOG(INFO) << "Moving average = " << mAvgTestInt0;
      if (isSMA) {
        // populate CCDB
      }
    }
  }
}

//______________________________________________________________________

void DCSProcessor::doSimpleMovingAverage(int nelements, std::deque<int>& vect, float& avg, bool& isSMA)
{

  // Do simple moving average on vector of ints
  if (vect.size() < nelements) {
    avg += vect[vect.size() - 1];
    return;
  }
  if (vect.size() == nelements) {
    avg += vect[vect.size() - 1];
    avg /= nelements;
    isSMA = true;
    return;
  }
  avg += (vect[vect.size() - 1] - vect[0]) / nelements;
  vect.pop_front();
  isSMA = true;
}

//______________________________________________________________________

void DCSProcessor::processDoubles()
{

  // function to process aliases of Double type
}

//______________________________________________________________________

void DCSProcessor::processUInts()
{

  // function to process aliases of UInt type
}

//______________________________________________________________________

void DCSProcessor::processBools()
{

  // function to process aliases of Bool type
}

//______________________________________________________________________

void DCSProcessor::processStrings()
{

  // function to process aliases of String type
}

//______________________________________________________________________

void DCSProcessor::processTimes()
{

  // function to process aliases of Time type
}

//______________________________________________________________________

void DCSProcessor::processBinaries()
{

  // function to process aliases of Time binary
}

//______________________________________________________________________

uint64_t DCSProcessor::processFlag(const uint64_t flags, const char* alias)
{

  // function to process the flag. the return code zero means that all is fine.
  // anything else means that there was an issue

  if (flags & DataPointValue::KEEP_ALIVE_FLAG) {
    LOG(INFO) << "KEEP_ALIVE_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::END_FLAG) {
    LOG(INFO) << "END_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::FBI_FLAG) {
    LOG(INFO) << "FBI_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::NEW_FLAG) {
    LOG(INFO) << "NEW_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::DIRTY_FLAG) {
    LOG(INFO) << "DIRTY_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::TURN_FLAG) {
    LOG(INFO) << "TURN_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::WRITE_FLAG) {
    LOG(INFO) << "WRITE_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::READ_FLAG) {
    LOG(INFO) << "READ_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::OVERWRITE_FLAG) {
    LOG(INFO) << "OVERWRITE_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::VICTIM_FLAG) {
    LOG(INFO) << "VICTIM_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::DIM_ERROR_FLAG) {
    LOG(INFO) << "DIM_ERROR_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::BAD_DPID_FLAG) {
    LOG(INFO) << "BAD_DPID_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::BAD_FLAGS_FLAG) {
    LOG(INFO) << "BAD_FLAGS_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::BAD_TIMESTAMP_FLAG) {
    LOG(INFO) << "BAD_TIMESTAMP_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::BAD_PAYLOAD_FLAG) {
    LOG(INFO) << "BAD_PAYLOAD_FLAG active for DP " << alias;
  }
  if (flags & DataPointValue::BAD_FBI_FLAG) {
    LOG(INFO) << "BAD_FBI_FLAG active for DP " << alias;
  }

  return 0;
}
