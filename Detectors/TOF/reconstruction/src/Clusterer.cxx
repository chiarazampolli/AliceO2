// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Clusterer.cxx
/// \brief Implementation of the TOF cluster finder
#include <algorithm>
#include "FairLogger.h" // for LOG
#include "DataFormatsTOF/Cluster.h"
#include "TOFReconstruction/Clusterer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::tof;

//__________________________________________________
Clusterer::Clusterer() 
{

  // empty for now
  
}
//__________________________________________________
void Clusterer::process(DataReader& reader, std::vector<Cluster>& clusters)
{
  reader.init();

  while (reader.getNextStripData(mStripData)) {
    LOG(DEBUG) << "TOFClusterer got Strip " << mStripData.stripID << " Ndigits "
               << mStripData.digits.size() << FairLogger::endl;
    
    processStrip(clusters);
  }
}

//__________________________________________________
void Clusterer::processStrip(std::vector<Cluster>& clusters)
{
  // method to clusterize the current strip

  Int_t detId[5];
  Int_t chan, chan2, chan3;
  Int_t strip1,strip2;
  Int_t iphi,iphi2,iphi3;
  Int_t ieta,ieta2,ieta3; // it is the number of padz-row increasing along the various strips
  
  for (int idig = 1; idig < mStripData.digits.size()-1; idig++) {
    Digit* dig = &mStripData.digits[idig];
    if (dig->isUsedInCluster()) continue; // the digit was already used to build a cluster

    dig->getPhiAndEtaIndex(iphi, ieta);

    // first we make a cluster out of the digit
    clusters.emplace_back();
    Cluster& c = clusters[noc];
    c.addContributingDigit(dig);
    
    for (int idigNext = idig+1; idigNext < mStripData.digits.size(); idigNext++) {
      Digit* digNext = &mStripData.digits[idigNext];
      if (digNext->isUsedInCluster()) continue; // the digit was already used to build a cluster
      // check if the TOF time are close enough to be merged; if not, it means that nothing else will contribute to the cluster (since digits are ordered in time)
      float timeDig = c.GetTime();
      float timeDigNext = digNext->GetTDC()*Geo::TDCBIN; // we assume it calibrated (for now); in ps
      if(timeDigNext - timeDig > 500/*in ps*/) break;
      digNext->getPhiAndEtaIndex(iphi2, ieta2);

      // check if the fired pad are close in space
      if((TMath::Abs(iphi-iphi2)>1) || (TMath::Abs(ieta-ieta2)>1)) 
	continue;

      // if we are here, the digit contributes to the cluster
      c.addContributingDigit(digNext);
      
    } // loop on the second digit

    c.buildCluster();

  } // loop on the first digit
}

 

