// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// example to run: o2-dpl-dcs-proxy --dcs-proxy '--channel-config "name=dcs-proxy,type=pull,method=connect,address=tcp://10.11.28.22:60000,rateLogging=1,transport=zeromq"' -b

#include "Framework/WorkflowSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/Lifetime.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/ExternalFairMQDeviceProxy.h"
#include <vector>

using namespace o2::framework;

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
  Outputs dcsOutput;
  dcsOutput.emplace_back("DCS", "DATAPOINTS", 0, Lifetime::Timeframe);
  
  DataProcessorSpec dcsProxy = specifyExternalFairMQDeviceProxy(
    "dcs-proxy",
    std::move(dcsOutput),
    "type=pull,method=connect,address=tcp://aldcsadaposactor:60000,rateLogging=1,transport=zmq",
    incrementalConverter(dcsOutput[0], 0, 1));

  WorkflowSpec workflow;
  workflow.emplace_back(dcsProxy);
  return workflow;
}
