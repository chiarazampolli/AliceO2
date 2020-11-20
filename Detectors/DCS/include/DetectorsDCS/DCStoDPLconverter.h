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
    
    static std::unordered_map<std::string, FairMQParts> outParts;
    static std::unordered_map<std::string, vector<DPCOM>> outputs;
    static std::unordered_map<DPID, DPCOM> cache;  // will keep only the latest measurement in the 1-second wide window for each DPID
    static auto timer = std::chrono::high_resolution_clock::now(); 

    // do we need to reserve the size of each element of outputs, according to how many DPs we expect at maximum per group?
    
    // We first iterate over the parts of the received message
    for (size_t i = 0; i < parts.Size(); ++i) { // DCS sends only 1 part, but we should be able to receive more
      std::unordered_map<std::string, FairMQParts> outParts; // if we want to add each group as a FairMQPart

      auto nDPCOM = (parts.At(i)->GetSize()) / sizeof(DPCOM); // number of DPCOM in current part
      for (size_t j = 0; j < nDPCOM; j++) { 
	const auto* src = reinterpret_cast<DPCOM*>(parts.At(i)->GetData() + id * sizeof(DPCOM));
	auto dst = cache[src->id]; // this is needed in case in the 1s window we get a new value for the same DP
	memcpy(&dst, src, sizeof(DPCOM));

	LOG(DEBUG) << "Reading from DCS: i = " << i << ", DPCOM = " << dst;

      }
    }

    auto timerNow = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<1>> duration = timerNow - timer;
    if (duration.count()>1) { //did we accumulate for 1 sec?
      *timesliceId += step; // we increment only if we send something     
      // do we need to set the time?
      // resetting timerNow to the moment when we send the output
	timerNow = timer;
      	for (auto& it : cache) {
	  // in the cache we have the final values of the DPs that we should put in the output
	  DPID dpid = it.first;
	  auto group = dpid2group[dpid]; // group to which this DP belongs
	  auto& groupOutput = outputs[group];
	  groupOutput.push_back(it.second);
	}
	for (auto& it : outputs) {
	  o2h::DataHeader hdr(it.first, "DCS", 0);
	  hdr.payloadSerializationMethod = o2h::gSerializationMethodNone;
	  hdr.splitPayloadParts = 1;
	  hdr.splitPayloadIndex = 1;
	  o2::header::Stack headerStack{hdr, o2::framework::DataProcessingHeader{*timesliceId, 0}};
	  memcpy(hdMessage->GetData(), headerStack.data(), headerStack.size());

	}
      }
    }


		 
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

