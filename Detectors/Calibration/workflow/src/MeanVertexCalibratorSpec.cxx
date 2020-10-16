// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   MeanVertexCalibratorSpec.cxx

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Logger.h"
#include "CalibrationWorkflow/MeanVertexCalibratorSpec.h"

using namespace o2::framework;

namespace o2
{
namespace calibration
{
void MeanVertexCalibDevice::init(InitContext& ic) {
 
}

//_____________________________________________________________

void MeanVertexCalibDevice::run(o2::framework::ProcessingContext& pc) {
}

//_____________________________________________________________

void MeanVertexCalibDevice::endOfStream(o2::framework::EndOfStreamContext& ec) {
}

//_____________________________________________________________

void MeanVertexCalibDevice::sendOutput(DataAllocator& output) {
}

//_____________________________________________________________

DataProcessorSpec getMeanVertexCalibDeviceSpec() {

  std::vector<OutputSpec> outputs;
  outputs.emplace_back(ConcreteDataTypeMatcher{clbUtils::gDataOriginCLB, clbUtils::gDataDescriptionCLBPayload});
  outputs.emplace_back(ConcreteDataTypeMatcher{clbUtils::gDataOriginCLB, clbUtils::gDataDescriptionCLBInfo});

  return DataProcessorSpec{
    "mean-vertex-calibration",
    outputs,
    Inputs{{"input", "GLO", "PVTX"}},
    AlgorithmSpec{adaptFromTask<MeanVertexCalibDevice>()},
      Options{}};

}

} // namespace calibration
} // namespace o2
