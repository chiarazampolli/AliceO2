// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TOFWorkflow/TOFClusterizerSpec.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Framework/Lifetime.h"
#include "Framework/Task.h"
#include "Headers/DataHeader.h"
#include "TOFReconstruction/Clusterer.h"
#include "TOFReconstruction/DataReader.h"
#include "DataFormatsTOF/Cluster.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "DataFormatsTOF/CalibLHCphaseTOF.h"
#include "DataFormatsTOF/CalibTimeSlewingParamTOF.h"
#include "TOFCalibration/CalibTOFapi.h"

#include <memory> // for make_shared, make_unique, unique_ptr
#include <vector>

using namespace o2::framework;

namespace o2
{
namespace tof
{

// use the tasking system of DPL
// just need to implement 2 special methods init + run (there is no need to inherit from anything)
class TOFDPLClustererTask
{
  using MCLabelContainer = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
  bool mUseMC = true;
  bool mUseCCDB = false;

 public:
  explicit TOFDPLClustererTask(bool useMC, bool useCCDB) : mUseMC(useMC), mUseCCDB(useCCDB) {}
  void init(framework::InitContext& ic)
  {
    // nothing special to be set up
  }

  void run(framework::ProcessingContext& pc)
  {
    static bool finished = false;
    if (finished) {
      return;
    }
    // get digit data
    auto digits = pc.inputs().get<std::vector<std::vector<o2::tof::Digit>>*>("tofdigits");
    auto labelvector = std::make_shared<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>>();
    if (mUseMC) {
      auto digitlabels = pc.inputs().get<std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>*>("tofdigitlabels");
      *labelvector.get() = std::move(*digitlabels);
      mClusterer.setMCTruthContainer(&mClsLabels);
      mClsLabels.clear();
    }

    o2::dataformats::CalibLHCphaseTOF lhcPhaseObj;
    o2::dataformats::CalibTimeSlewingParamTOF channelCalibObj;

    if(mUseCCDB){ // read calibration objects from ccdb
      // check LHC phase
      auto lhcPhase = pc.inputs().get<o2::dataformats::CalibLHCphaseTOF*>("tofccdbLHCphase");
      printf("\n\n\n\n\n\n\n\nCCDB: lhcPhase size = %d\n\n\n\n\n\n\n\n", lhcPhase->size());
      auto channelCalib = pc.inputs().get<o2::dataformats::CalibTimeSlewingParamTOF*>("tofccdbChannelCalib");
      printf("\n\n\n\n\n\n\n\nCCDB: channelCalib size = %d\n\n\n\n\n\n\n\n", channelCalib->size());

      o2::dataformats::CalibLHCphaseTOF lhcPhaseObjTmp = std::move(*lhcPhase);
      o2::dataformats::CalibTimeSlewingParamTOF channelCalibObjTmp = std::move(*channelCalib);


      // make a copy in global scope
      lhcPhaseObj = lhcPhaseObjTmp;
      channelCalibObj = channelCalibObjTmp;
    }
    else{ // calibration objects set to zero
      lhcPhaseObj.addLHCphase(0, 0);
      lhcPhaseObj.addLHCphase(2000000000, 0);

      for (int ich = 0; ich <o2::dataformats::CalibTimeSlewingParamTOF::NCHANNELS; ich++){
	  channelCalibObj.addTimeSlewingInfo(ich, 0, 0);
	  int sector = ich/o2::dataformats::CalibTimeSlewingParamTOF::NCHANNELXSECTOR;
	  int channelInSector = ich%o2::dataformats::CalibTimeSlewingParamTOF::NCHANNELXSECTOR;
	  channelCalibObj.setFractionUnderPeak(sector, channelInSector, 1);
      }
    }

    printf("\n\n\n\n\n\n\n\nlhcPhaseObj size = %d\n\n\n\n\n\n\n\n", lhcPhaseObj.size());
    printf("\n\n\n\n\n\n\n\nchannelCalibObj size = %d\n\n\n\n\n\n\n\n", channelCalibObj.size());

    o2::tof::CalibTOFapi calibapi(long(0), &lhcPhaseObj, &channelCalibObj);

    mClusterer.setCalibApi(&calibapi);

    // call actual clustering routine
    mClustersArray.clear();

    for (int i = 0; i < digits->size(); i++) {
      printf("# TOF readout window for clusterization = %i\n", i);
      auto digitsRO = digits->at(i);
      mReader.setDigitArray(&digitsRO);
      if (mUseMC) {
        mClusterer.process(mReader, mClustersArray, &(labelvector->at(i)));
      } else
        mClusterer.process(mReader, mClustersArray, nullptr);
    }
    LOG(INFO) << "TOF CLUSTERER : TRANSFORMED " << digits->size()
              << " DIGITS TO " << mClustersArray.size() << " CLUSTERS";

    // send clusters
    pc.outputs().snapshot(Output{"TOF", "CLUSTERS", 0, Lifetime::Timeframe}, mClustersArray);
    // send labels
    if (mUseMC)
      pc.outputs().snapshot(Output{"TOF", "CLUSTERSMCTR", 0, Lifetime::Timeframe}, mClsLabels);

    // declare done
    finished = true;
    //pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    pc.services().get<ControlService>().endOfStream();
  }

 private:
  DigitDataReader mReader; ///< Digit reader
  Clusterer mClusterer;    ///< Cluster finder

  std::vector<Cluster> mClustersArray; ///< Array of clusters
  MCLabelContainer mClsLabels;
};

o2::framework::DataProcessorSpec getTOFClusterizerSpec(bool useMC, bool useCCDB)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("tofdigits", "TOF", "DIGITS", 0, Lifetime::Timeframe);
  if(useCCDB){
    inputs.emplace_back("tofccdbLHCphase", "TOF", "LHCphase");
    inputs.emplace_back("tofccdbChannelCalib", "TOF", "ChannelCalib");
  }
  if (useMC)
    inputs.emplace_back("tofdigitlabels", "TOF", "DIGITSMCTR", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    "TOFClusterer",
    inputs,
    Outputs{OutputSpec{"TOF", "CLUSTERS", 0, Lifetime::Timeframe},
            OutputSpec{"TOF", "CLUSTERSMCTR", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<TOFDPLClustererTask>(useMC, useCCDB)},
    Options{/* for the moment no options */}};
}

} // end namespace tof
} // end namespace o2
