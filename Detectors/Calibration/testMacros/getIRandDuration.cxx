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

// executable to get the interaction rate and duration of a run from CCDB

#include <fstream>
#include <stdio.h>
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Scalers.h"
#include "DataFormatsCTP/Configuration.h"
#include <boost/program_options.hpp>
#endif

using namespace o2::ctp;
namespace bpo = boost::program_options;

const double orbitDuration = 88.924596234; // us

void writeIRtoFile(float ir)
{

  FILE* fptr = fopen("IR.txt", "w");
  if (fptr == NULL) {
    printf("ERROR: Could not open file to write IR!");
    return;
  }
  fprintf(fptr, "%.2f", ir);
  fclose(fptr);
}

void writeDurationToFile(long duration)
{

  FILE* fptr = fopen("Duration.txt", "w");
  if (fptr == NULL) {
    printf("ERROR: Could not open file to write IR!");
    return;
  }
  fprintf(fptr, "%ld", duration);
  fclose(fptr);
}

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[], bpo::variables_map& vm)
{
  options.add_options()(
    "run,r", bpo::value<int>()->default_value(8), "Run number to inspect")(
    "enable-debug,d", bpo::value<bool>()->default_value(false)->implicit_value(true), "Enable debug logs")(
    "help,h", "Produce help message.");

  try {
    bpo::store(parse_command_line(argc, argv, options), vm);

    // help
    if (vm.count("help")) {
      std::cout << options << std::endl;
      return false;
    }

    bpo::notify(vm);
  } catch (const bpo::error& e) {
    std::cerr << e.what() << "\n\n";
    std::cerr << "Error parsing command line arguments; Available options:\n";

    std::cerr << options << std::endl;
    return false;
  }
  return true;
}

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  bpo::variables_map vm;
  if (!initOptionsAndParse(options, argc, argv, vm)) {
    return -1;
  }

  auto run = vm["run"].as<int>();
  auto debug = vm["enable-debug"].as<bool>();

  float ir = 0.f;
  long duration = 0;
  // duration as O2end - O2start:
  auto& ccdb_inst = o2::ccdb::BasicCCDBManager::instance();
  ccdb_inst.setURL("https://alice-ccdb.cern.ch");
  std::pair<uint64_t, uint64_t> run_times = ccdb_inst.getRunDuration(run);
  long run_duration = long(run_times.second - run_times.first);

  LOGP(info, "Checking IR anbd duration");
  if (run < 523141) {
    // LHC22c, d, e, f
    LOGP(info, "Run number < 523141 --> we are in 22c, d, e, or f, so IR is < 100 kHz, writing 0.f");
    LOGP(info, "In addition, the duration for these runs is O2end - O2start: if the run was short, this might overestimate the duration");
    // In these runs, sometimes the CCDB does not contain correct scalers, so we use 0 as a placeholder
    writeIRtoFile(ir);
    writeDurationToFile(duration);
    return 0;
  }

  o2::ccdb::CcdbApi ccdb_api;
  ccdb_api.init("https://alice-ccdb.cern.ch");
  // access SOR and EOR timestamps
  std::map<string, string> headers, metadata;
  headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadata, -1);
  int64_t tsSOR = atol(headers["SOR"].c_str()); // ms
  int64_t tsEOR = atol(headers["EOR"].c_str()); // ms
  LOGP(info, "tsSOR = {} ms, tsEOR = {} ms", tsSOR, tsEOR);

  // Extract CTP info
  std::map<std::string, std::string> metadataCTP;
  metadataCTP["runNumber"] = Form("%d", run);
  ccdb_inst.setFatalWhenNull(false);
  o2::ctp::CTPRunScalers* scl = ccdb_inst.getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", tsSOR, metadataCTP);
  if (!scl) {
    LOGP(info, "CTP/Calib/Scalers object does not exist in production CCDB, trying test CCDB");
    ccdb_inst.setURL("http://ccdb-test.cern.ch:8080");
    scl = ccdb_inst.getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", tsSOR, metadata);
    if (!scl) {
      LOGP(info, "Cannot get IR for run {} neither from production nor test CCDB, writing -1.f", run);
      LOGP(info, "In addition, the duration for these runs is O2end - O2start: if the run was short, this might overestimate the duration");
      ir = -1.f;
      writeIRtoFile(ir);
      writeDurationToFile(duration);
      return 0;
    }
  }

  scl->convertRawToO2();
  std::vector<CTPScalerRecordO2> mScalerRecordO2 = scl->getScalerRecordO2();
  int n = mScalerRecordO2.size();
  if (n != 0) {
    std::int64_t totScalers = 0;
    std::vector<int64_t> vOrbit;
    std::vector<int64_t> vScaler;
    int i = 0;
    for (auto& record : mScalerRecordO2) {
      if (debug) {
        record.printStream(std::cout);
      }
      std::vector<CTPScalerO2>& scalers = record.scalers;
      o2::InteractionRecord& intRecord = record.intRecord;
      vOrbit.push_back(intRecord.orbit);
      if (debug) {
        LOGP(info, "{} orbit = {} scalers = {}", i, intRecord.orbit, scalers[0].lmBefore);
      }
      vScaler.push_back(scalers[0].lmBefore); // use scalers for class 0 (usually TVX). TODO: extract info on class id from trigger config
      totScalers += scalers[0].lmBefore;
      ++i;
    }

    duration = std::round((vOrbit.back() - vOrbit.front()) * orbitDuration * 1e-6); // s
    ir = float(vScaler.back() - vScaler.front()) / duration;
    LOGP(info, "run {}: orbit.front = {} orbit.back = {} duration = {} s scalers = {} IR = {} Hz", run, vOrbit.front(), vOrbit.back(), duration, vScaler.back() - vScaler.front(), ir);
  }

  if (ir < 100000) {
    LOGP(info, "IR < 100 kHz");
  } else {
    LOGP(info, "IR > 100 kHz");
  }
  writeIRtoFile(ir);
  writeDurationToFile(duration);
  return 0;
}
