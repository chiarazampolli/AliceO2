#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <sstream>

#include <TStopwatch.h>

#include "FairLogger.h"
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairSystemInfo.h"

#include "TOFSimulation/DigitizerTask.h"
#endif

void run_digi_tof(Int_t nEvents = 10, Float_t rate=50.e5)
{
  // if rate>0 then continuous simulation for this rate will be performed

  // Initialize logger
  // FairLogger* logger = FairLogger::GetLogger();
  // logger->SetLogVerbosityLevel("LOW");
  // logger->SetLogScreenLevel("DEBUG");

  // Input and output file name
  std::stringstream inputfile, outputfile, paramfile;
  //inputfile << "AliceO2_" << mcEngine << ".mc_" << nEvents << "_event.root";
  //paramfile << "AliceO2_" << mcEngine << ".params_" << nEvents << ".root";
  //outputfile << "AliceO2_" << mcEngine << ".digi_" << nEvents << "_event.root";
  inputfile << "o2sim.root";
  paramfile << "o2sim_par.root";
  outputfile << "o2sim_digi.root";

  // Setup timer
  TStopwatch timer;

  // Setup FairRoot analysis manager
  FairRunAna* fRun = new FairRunAna();
  FairFileSource* fFileSource = new FairFileSource(inputfile.str().c_str());
  fRun->SetSource(fFileSource);
  fRun->SetOutputFile(outputfile.str().c_str());

  if (rate > 0) {
    fFileSource->SetEventMeanTime(1.e9 / rate); // is in us
  }

  // Setup Runtime DB
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(paramfile.str().c_str());
  rtdb->setFirstInput(parInput1);

  // Setup digitizer
  o2::tof::DigitizerTask* digi = new o2::tof::DigitizerTask();
  digi->setContinuous(rate > 0);
//  digi->setFairTimeUnitInNS(1.0); // tell in which units (wrt nanosecond) FAIT timestamps are
  fRun->AddTask(digi);

  fRun->Init();

  timer.Start();
  fRun->Run();

}
