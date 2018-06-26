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

    mNumberOfContributingDigits = 0;
    dig->getPhiAndEtaIndex(iphi, ieta);

    // first we make a cluster out of the digit
    clusters.emplace_back();
    Cluster& c = clusters[noc];
    addContributingDigit(dig);
    
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
      addContributingDigit(digNext);
      
    } // loop on the second digit

    buildCluster(c);

  } // loop on the first digit
}
//______________________________________________________________________
void Cluster::addContributingDigit(Digit* dig) {

  // adding a digit to the array that stores the contributing ones

  if (mNumberOfContributingDigits == 6) {
    LOG(ERROR) << "The cluster has already 6 digits associated to it, we cannot add more; returning without doing anything" << FairLogger::endl;
  }
  mContributingDigit[mNumberOfContributingDigits] = dig;
  mNumberOfContributingDigits++;
  dig->setIsUsedInCluster();

  return;
}

//_____________________________________________________________________
void Cluster::buildCluster(Cluster& c) {

  // here we finally build the cluster from all the digits contributing to it

  Digit* temp;
  for (idig = 1; idig < mNumberOfContributingDigits; idig++){
    // the digit[0] will be the main one
    if (mContributingDigit[idig]->getTOT() > mContributingDigit[0]->getTOT()){
      temp = mContributingDigit[0];
      mContributingDigit[0] = mContributingDigit[idig];
      mContributingDigit[idig] = temp;
    }
  }

  c.setMainContributingChannel(mContributingDigit[0]->getChannel());
  c.setTime(mContributingDigit[0]->getTDC()*Geo::TDCBIN); // time in ps (for now we assume it calibrated)
  c.setTot(mContributingDigit[0]->getTOT()*Geo::TOTBIN*1E-3); // TOT in ns (for now we assume it calibrated)
  //setL0L1Latency(); // to be filled (maybe)
  //setDeltaBC(); // to be filled (maybe)

  int chan1, chan2;
  int phi1, phi2;
  int eta1, eta2;
  int deltaPhi, deltaEta;
  int mask;

  mContributingDigit[0]->getPhiAndEtaIndex(phi1, eta1);
  // now set the mask with the secondary digits
  for (idig = 1; idig < mNumberOfContributingDigits; idig++){
    mContributingDigit[idig]->getPhiAndEtaIndex(phi2, eta2);
    deltaPhi = phi1-phi2;
    deltaEta = eta1-eta2;
    if (deltaPhi == 1) { // the digit is to the LEFT of the cluster; let's check about UP/DOWN/Same Line
      if (deltaEta == 1) { // the digit is DOWN LEFT wrt the cluster
	mask = kDownLeft;
      }
      else if (deltaEta == -1) { // the digit is UP LEFT wrt the cluster
	mask = kUpLeft;
      }
      else { // the digit is LEFT wrt the cluster
	mask = kLeft;
      }
    }
    else if (deltaPhi == -1) { // the digit is to the RIGHT of the cluster; let's check about UP/DOWN/Same Line
      if (deltaEta == 1) { // the digit is DOWN RIGHT wrt the cluster
	mask = kDownRight;
      }
      else if (deltaEta == -1) { // the digit is UP RIGHT wrt the cluster
	mask = kUpRight;
      }
      else { // the digit is RIGHT wrt the cluster
	mask = kRight;
      }
    }
    else if (deltaPhi == 0) { // the digit is on the same column as the cluster; is it UP or Down?
      if (deltaEta == 1) { // the digit is DOWN wrt the cluster
	mask = kDown;
      }
      else if (deltaEta == -1) { // the digit is UP wrt the cluster
	mask = kUp;
      }
      else { // impossible!!
	LOG(DEBUG) << " Check what is going on, the digit you are trying to merge to the cluster must be in a different channels... " << FairLogger::endl;
      }
    }
    else { // impossible!!! We checked above... 
      LOG(DEBUG) << " Check what is going on, the digit you are trying to merge to the cluster is too far from the cluster, you should have not got here... " << FairLogger::endl;
    }
    c.addBitInContributingChannels(mask);
  }

  return;
  
}



