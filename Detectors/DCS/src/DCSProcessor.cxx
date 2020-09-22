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

using namespace o2::dcs;

using ints = std::vector<int>;
using chars = std::vector<char>;
using doubles = std::vector<double>;
using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DPCOM = o2::dcs::DataPointCompositeObject;
using DeliveryType = o2::dcs::DeliveryType;

//ClassImp(o2::dcs::DCSProcessor);

int DCSProcessor::process(std::unordered_map<DPID, DPVAL> map) {

  // process function to do "something" with the DCS map that is passed

  // first, we need to check if there are the Data Points that we need

  //for (auto it = std::begin(mAliaseschars); it != std::end(mAliaseschars); ++it) {

  for (auto i = 0; i != mAliaseschars.size(); ++i) {
    if (map.find(mAliaseschars[i]) == map.end()) {
      LOG(ERROR) << "Element " << mAliaseschars[i] << " not found " << std::endl;
      return 1;
    }
    DeliveryType tt = mAliaseschars[i].get_type();
    if (tt != DeliveryType::RAW_CHAR) {
      LOG(FATAL) << "Delivery Type of alias " << mAliaseschars[i].get_alias() << " does not match definition in DCSProcessor! Please fix";
    }    
    DPVAL val = map[mAliaseschars[i]];
    mDpschars[i].push_back(val.payload_pt1);
  }

  for (auto i = 0; i != mAliasesints.size(); ++i) {
    if (map.find(mAliasesints[i]) == map.end()) {
      LOG(ERROR) << "Element " << mAliasesints[i] << " not found " << std::endl;
      return 1;
    }
    DeliveryType tt = mAliasesints[i].get_type();
    if (tt != DeliveryType::RAW_INT) {
      LOG(FATAL) << "Delivery Type of alias " << mAliasesints[i].get_alias() << " does not match definition in DCSProcessor! Please fix";
    }    
    DPVAL val = map[mAliasesints[i]];
    mDpsints[i].push_back(val.payload_pt1);
  }

  for (auto i = 0; i != mAliasesdoubles.size(); ++i) {
    if (map.find(mAliasesdoubles[i]) == map.end()) {
      LOG(ERROR) << "Element " << mAliasesdoubles[i] << " not found " << std::endl;
      return 1;
    }
    DeliveryType tt = mAliasesdoubles[i].get_type();
    if (tt != DeliveryType::RAW_DOUBLE) {
      LOG(FATAL) << "Delivery Type of alias " << mAliasesdoubles[i].get_alias() << " does not match definition in DCSProcessor! Please fix";
    }    
    DPVAL val = map[mAliasesdoubles[i]];
    mDpsdoubles[i].push_back(val.payload_pt1);
  }

    /*
    DPVAL val = map[*it];
    
    if (tt == DeliveryType::RAW_INT) {
      mDps[*it].push_back(std::get<int>(val.payload_pt1));
    }
    else if (tt == DeliveryType::RAW_DOUBLE) {
      mDps[*it].push_back(std::get<double>(val.payload_pt1));
      //mDps[*it].push_back((double)val.payload_pt1);
    }
    else if (tt == DeliveryType::RAW_CHAR) {
      mDps[*it].push_back(std::get<char>(val.payload_pt1));
      //mDps[*it].push_back((char)val.payload_pt1);
    } 
  
  }
    */

  return 0;
}
