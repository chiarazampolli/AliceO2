// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// example to run:
// o2-dcs-dcs-proxy --dcs-proxy '--channel-config "name=dcs-proxy,type=pull,method=connect,address=tcp://10.11.28.22:60000,rateLogging=1,transport=zeromq"' -b

#include "Framework/WorkflowSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/Lifetime.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/ExternalFairMQDeviceProxy.h"
#include "DetectorsDCS/DataPointIdentifier.h"
#include "DetectorsDCS/DataPointValue.h"
#include "DetectorsDCS/DeliveryType.h"
#include "DCStoDPLconverter.h"
#include <vector>
#include <unordered_map>

using namespace o2::framework;
using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DeliveryType = o2::dcs::DeliveryType;

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  /*  workflowOptions.push_back(
    ConfigParamSpec{
      "dataspec", VariantType::String, "A:FLP/RAWDATA;B:FLP/DISTSUBTIMEFRAME/0", {"selection string for the data to be proxied"}});

  workflowOptions.push_back(
    ConfigParamSpec{
      "throwOnUnmatched", VariantType::Bool, false, {"throw if unmatched input data is found"}});
  */
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{

  /*
  std::unordered_map<DPID, OutputSpec> dp2OutputSpec;
  OutputSpec osTestDouble(ConcreteDataTypeMatcher{"DCS", "TESTDOUBLES"});
  OutputSpec osTestBool(ConcreteDataTypeMatcher{"DCS", "TESTBOOLS"});

  // this should be filled e.g. reading from CCDB. now for tests I hardcode it
  DPID dpidtmp;
  std::string dpAlias = "ADAPOS_LG/TEST_000260";
  DPID::FILL(dpidtmp, dpAlias, DeliveryType::RAW_DOUBLE);
  dp2OutputSpec[dpidtmp] = osTestDouble;

  dpAlias = "ADAPOS_LG/TEST_000261";
  DPID::FILL(dpidtmp, dpAlias, DeliveryType::RAW_DOUBLE);
  dp2OutputSpec[dpidtmp] = osTestDouble;

  dpAlias = "ADAPOS_LG/TEST_000274";
  DPID::FILL(dpidtmp, dpAlias, DeliveryType::RAW_BOOL);
  dp2OutputSpec[dpidtmp] = osTestBool;

  dpAlias = "ADAPOS_LG/TEST_000275";
  DPID::FILL(dpidtmp, dpAlias, DeliveryType::RAW_BOOL);
  dp2OutputSpec[dpidtmp] = osTestBool;
  */

  DPID dpidtmp;

  std::unordered_map<DPID, o2h::DataDescription> dpid2DataDesc;
  std::string dpAlias = "ADAPOS_LG/TEST_000261";
  DPID::FILL(dpidtmp, dpAlias, DeliveryType::RAW_DOUBLE);
  dpid2DataDesc[dpidtmp] = "COMMON"; // i.e. this will go to {DCS/COMMON/0} OutputSpec

  // RS: here we should complete the attribution of different DPs to different outputs
  // ...

  // now collect all required outputs to define OutputSpecs for specifyExternalFairMQDeviceProxy
  std::unordered_map<o2h::DataDescription, int, std::hash<o2h::DataDescription>> outMap;
  for (auto itdp : dpid2DataDesc) {
    outMap[itdp.second]++;
  }

  Outputs dcsOutputs;
  for (auto itout : outMap) {
    dcsOutputs.emplace_back("DCS", itout.first, 0, Lifetime::Timeframe);
  }

  DataProcessorSpec dcsProxy = specifyExternalFairMQDeviceProxy(
    "dcs-proxy",
    std::move(dcsOutputs),
    "type=pull,method=connect,address=tcp://aldcsadaposactor:60000,rateLogging=1,transport=zeromq",
    dcs2dpl(dpid2DataDesc, 0, 1));

  WorkflowSpec workflow;
  workflow.emplace_back(dcsProxy);
  return workflow;
}
