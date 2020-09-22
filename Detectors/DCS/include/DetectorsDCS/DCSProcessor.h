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
#include <variant>
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
  
  using ints = std::vector<int>;
  using chars = std::vector<char>;
  using doubles = std::vector<double>;
  using DPID = o2::dcs::DataPointIdentifier;
  using DPVAL = o2::dcs::DataPointValue;
  using DPCOM = o2::dcs::DataPointCompositeObject;
 
public:
  DCSProcessor() = default;
  ~DCSProcessor() = default;

  void init(std::vector<DPID> aliaseschars, std::vector<DPID> aliasesints, std::vector<DPID> aliasesdoubles) {

    for (auto it = std::begin(aliaseschars); it != std::end(aliaseschars); ++it){
      //mAliaseschars.emplace_back(*it->get_alias(), *it->get_type());
      if ((*it).get_type() != DeliveryType::RAW_CHAR) {
	LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a char";
      }
      mAliaseschars.emplace_back((*it).get_alias(), DeliveryType::RAW_CHAR);
    }
    for (auto it = std::begin(aliasesints); it != std::end(aliasesints); ++it){
      if ((*it).get_type() != DeliveryType::RAW_INT) {
	LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a int";
      }
      mAliasesints.emplace_back((*it).get_alias(), DeliveryType::RAW_INT);
    }
    for (auto it = std::begin(aliasesdoubles); it != std::end(aliasesdoubles); ++it){
      if ((*it).get_type() != DeliveryType::RAW_DOUBLE) {
	LOG(FATAL) << "Type for data point " << *it << " does not match with expectations! It should be a double";
      }
      mAliasesdoubles.emplace_back((*it).get_alias(), DeliveryType::RAW_DOUBLE);
    }

    mDpschars.resize(mAliaseschars.size());
    mDpsints.resize(mAliasesints.size());
    mDpsdoubles.resize(mAliasesdoubles.size());

  }
  
  int process(std::unordered_map<DPID, DPVAL> map);

 private:
  std::vector<chars> mDpschars;
  std::vector<ints> mDpsints;
  std::vector<doubles> mDpsdoubles;
  //  std::unordered_map<DPID, std::vector<DPVAL>> mDps;
  std::vector<DPID> mAliaseschars;
  std::vector<DPID> mAliasesints;
  std::vector<DPID> mAliasesdoubles;
  
  ClassDefNV(DCSProcessor, 0);
};

} // namespace dcs
} // namespace o2

#endif
