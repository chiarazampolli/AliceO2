// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   Utils.h
/// @brief  Utils and constants for calibration and related workflows

#ifndef O2_CALIBRATION_CONVENTIONS_H
#define O2_CALIBRATION_CONVENTIONS_H

#include <typeinfo>
#include <utility>
#include <fstream>
#include <TMemFile.h>
#include "Headers/DataHeader.h"
#include "CommonUtils/StringUtils.h"
#include "CCDB/CCDBTimeStampUtils.h"
#include "CCDB/CcdbObjectInfo.h"
#include "CommonUtils/MemFileHelper.h"
#include "CCDB/CcdbApi.h"
#include <vector>

namespace o2
{
namespace calibration
{

struct Utils {

  enum ValueType { Invalid = 0,
                   Interpolation,
                   ClosestAvailableFromBelow,
                   ClosestAvailableFromAbove,
                   SameAsRequested };

  static constexpr o2::header::DataOrigin gDataOriginCDBPayload{"CLP"}; // generic DataOrigin for calibrations payload
  static constexpr o2::header::DataOrigin gDataOriginCDBWrapper{"CLW"}; // generic DataOrigin for calibrations wrapper
  template <typename T>
  static void prepareCCDBobjectInfo(T& obj, o2::ccdb::CcdbObjectInfo& info, const std::string& path,
                                    const std::map<std::string, std::string>& md, long start, long end = -1);

  static std::pair<int, double> findPair(std::vector<std::pair<uint64_t, double>> vect, uint64_t timestamp)
  {
    // function to find the pair in the vector with the timestamp closest in the past to the one passed as argument
    // if two entries in the vector have the same timestamp, the first that is found is used
    std::vector<uint64_t> times;
    for (const auto& el : vect) {
      times.push_back(el.first);
    }
    auto lower = std::lower_bound(times.begin(), times.end(), timestamp);
    if (*lower == timestamp) {
      LOG(info) << "We found the element for the exact timestamp";
      return std::make_pair(SameAsRequested, vect[std::distance(times.begin(), lower)].second);
    } else if (lower == times.end()) {
      LOG(info) << "All values are smaller than the queried one " << timestamp << ", we return the closest available from below: " << times.back();
      return std::make_pair(ClosestAvailableFromBelow, vect.back().second);
    } else if (lower == times.begin()) {
      LOG(info) << "All values are greater that the queried one " << timestamp << ", we return the closest available from above: " << *times.begin();
      return std::make_pair(ClosestAvailableFromAbove, vect.begin()->second);
    } else {
      // doing interpolation
      const auto& p1 = vect[std::distance(times.begin(), lower) - 1];
      const auto& p2 = vect[std::distance(times.begin(), lower)];
      auto t1 = p1.first;
      auto t2 = p2.first;
      auto val1 = p1.second;
      auto val2 = p2.second;
      if (t1 == t2) {
        LOG(info) << "times are the same, cannot interpolate, returning closest from below";
        return std::make_pair(ClosestAvailableFromBelow, val1);
      }
      // in case values are not ordered in timestamp
      double val = t2 > t1 ? (timestamp - t1) * (val2 - val1) / (t2 - t1) + val1 : (timestamp - t1) * (val1 - val2) / (t1 - t2) + val1;
      LOG(info) << "Doing interpolation between (" << t1 << ", " << val1 << ") and (" << t2 << ", " << val2 << ") for t = " << timestamp << " --> " << val;
      return std::make_pair(Interpolation, val);
    }
    LOG(info) << "Something went wrong!";
    return std::make_pair(Invalid, 0);
  }
};

template <typename T>
void Utils::prepareCCDBobjectInfo(T& obj, o2::ccdb::CcdbObjectInfo& info, const std::string& path,
                                  const std::map<std::string, std::string>& md, long start, long end)
{

  // prepare all info to be sent to CCDB for object obj
  auto clName = o2::utils::MemFileHelper::getClassName(obj);
  auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
  info.setPath(path);
  info.setObjectType(clName);
  info.setFileName(flName);
  info.setStartValidityTimestamp(start);
  info.setEndValidityTimestamp(end);
  info.setMetaData(md);
}

} // namespace calibration
} // namespace o2

#endif
