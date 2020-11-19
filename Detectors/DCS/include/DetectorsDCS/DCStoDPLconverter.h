// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#ifndef O2_DCS_TO_DPL_CONVERTER
#define O2_DCS_TO_DPL_CONVERTER

#include "Framework/DataSpecUtils.h"
#include "Framework/ExternalFairMQDeviceProxy.h"
#include <fairmq/FairMQParts.h>
#include <fairmq/FairMQDevice.h>
#include "DetectorsDCS/DataPointIdentifier.h"
#include "DetectorsDCS/DataPointValue.h"
#include <unordered_map>

namespace o2
{
namespace dcs
{
using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DPCOM = o2::dcs::DataPointCompositeObject;
  

InjectorFunction dcs2dpl(std::vector<OutputSpec> specs, std::unordered_map<DPID, std::string>& dpid2group,
			 uint64_t startTime, uint64_t step) 
{

  auto timesliceId = std::make_shared<size_t>(startTime);

  return [specs, dpid2group, timesliceId, step](FairMQDevice& device, FairMQParts& parts, ChannelRetriever channelRetriever) { // why do we capture by copy?
    
    // We first iterate over the parts of the received message
    for (size_t i = 0; i < parts.Size(); ++i) {
      std::unordered_map<std::string, FairMQParts> outPart; // if we want to add each DP as a FairMQPart
      // now I have to build the message that goes to each OutputSpec
      // I need to select what I am interested in from the DCS data
      // I look in the DCS data, and I associate each DP to the message for the corresponding group
      DPCOM dpcom;
      auto nDPCOM = (parts.At(i)->GetSize()) / sizeof(DPCOM); // number of DPCOM in current part
      for (int i = 0; i < nDPCOM; i++) {
	memcpy(&dpcom, parts.At(i)->GetData() + i * sizeof(DPCOM), sizeof(DPCOM));
	LOG(DEBUG) << "Reading from DCS: i = " << i << ", DPCOM = " << dptmp;
	DPID dpid = dpcom.id;
	auto group = dpid2group[dpid]; // group to which this DP belongs
	// how to I create the message of the group, that will then go to the correct outputSpec?
	outPart[group].AddPart(std::move(...)); // ?  
	memcpy(...); // ?	  



		 
      for (spec : specs) { // loop over OutputSpec
	// create the data header of the message
	DataHeader dh;
	ConcreteDataMatcher matcher = DataSpecUtils::asConcreteDataMatcher(spec);
	dh.dataOrigin = matcher.origin;
	dh.dataDescription = matcher.description;
	dh.subSpecification = matcher.subSpec;	auto nDPs = ndpsPerOS[spec]; // how many DPs are associated to the OutputSpec	
	dh.payloadSize = nDPs * (sizeof(DPID) + sizeof(DPVAL)); // FIXME: this is wrong! only after we read the data that go to the OutputSpec we know the size --> should be moved later

	// data processing header 
	DataProcessingHeader dph{*timesliceId, 0};
	*timesliceId += step;
	
	//we have to move the incoming data
	o2::header::Stack headerStack{dh, dph};

	  
	// parts.At(i) is what DCS has sent me. Now I need to extract the data for the DPs corresponding to the OutputSpec ospairs[iospec].first
	sendOnChannel(device, std::move(headerStack), std::move(specPart), spec, channelRetriever);
      }
    }
      
  }
}
	     
} // namespace dcs
} // namespace o2

#endif /* O2_DCS_TO_DPL_CONVERTER_H */

