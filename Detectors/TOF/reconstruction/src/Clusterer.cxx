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
Clusterer::Clusterer() : mCurr(mColumn2 + 1), mPrev(mColumn1 + 1)
{
  std::fill(std::begin(mColumn1), std::end(mColumn1), -1);
  std::fill(std::begin(mColumn2), std::end(mColumn2), -1);
}
//__________________________________________________
void Clusterer::process(DataReader& reader, std::vector<Cluster>& clusters)
{
  reader.init();

  while (reader.getNextStripData(mStripData)) {
    LOG(DEBUG) << "TOFClusterer got Strip " << mStripData.stripID << " Ndigits "
               << mStripData.digits.size() << FairLogger::endl;
    
    initStrip();
    /*    for (int idig = 1; idig < mStripData.digits.size(); idig++)
      updateStrip(idig);
    */
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
    chan = dig->getChannel(); // note that inside the strip the digits are ordered per channel number
    Geo::getVolumeIndices(chan, detId); // Get volume index from channel index
    ieta = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
    if(detId[1]/*module*/ == 0) ieta += 0;
    else if(detId[1] == 1) ieta += 38;
    else if(detId[1] == 2) ieta += 76;
    else if(detId[1] == 3) ieta += 106;
    else if(detId[1] == 4) ieta += 144;
    iphi = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;

    // first we make a cluster out of the digit
    clusters.emplace_back();
    Cluster& c = clusters[noc];
    c.setTime(dig->getTDC()*Geo::TDCBIN); // time in ps (for now we assume it calibrated)
    c.setTot(dig->getTOT()*Geo::TOTBIN*1E-3); // TOT in ns (for now we assume it calibrated)
    //c.setL0L1Latency(); // to be filled (maybe)
    //c.setDeltaBC(); // to be filled (maybe)
    c.setMainContributingChannel(chan1);
    
    for (int idigNext = idig+1; idigNext < mStripData.digits.size(); idigNext++) {
      Digit* digNext = &mStripData.digits[idigNext];
      if (digNext->isUsedInCluster()) continue; // the digit was already used to build a cluster
      chan2 = digNext->getChannel(); // note that inside the strip the digits are ordered per channel number
      if (chan == chan2) {
	// the digits belong to the same channel, but were not merged; let's assume that they can make a cluster
	LOG(DEBUG) << " found two digits from the same channel that were not merged - they can still make a cluster " << FairLogger::endl;
      }
      Geo::getVolumeIndices(chan2, detId); // Get volume index from channel index
      ieta2 = detId[2]/*strip*/*2 + detId[3]/*pad Z*/;
      if(detId[1]/*module*/ == 0) ieta2 += 0;
      else if(detId[1] == 1) ieta2 += 38;
      else if(detId[1] == 2) ieta2 += 76;
      else if(detId[1] == 3) ieta2 += 106;
      else if(detId[1] == 4) ieta2 += 144;
      iphi2 = detId[0]/*phi sector*/*48 + detId[4]/*pad x*/;

      // if(ieta3-ieta > 2) j = fN; // because cluster are order along Z so you can skip all the rest, go to the next one ("i+1")

      // check if the fired pad are close in space
      if((TMath::Abs(iphi-iphi2)>1) || (TMath::Abs(ieta-ieta2)>1)) 
	continue;
      
      // check if the TOF time are close enough to be merged
      float timeDig = c.GetTime();
      float timeDigNext = digNext->GetTDC()*Geo::TDCBIN; // we assume it calibrated (for now); in ps
      
      if(TMath::Abs(timeDig - timeDigNext) > 500/*in ps*/) continue;

      int deltaPhi = iphi - iphi2; // difference in x between cluster and digit to be merged
      int deltaEta = ieta - ieta2; // difference in z between cluster and digit to be merged
      // merge the digit to the existing cluster
      bool swapped = mergeClusters(c, digNext, deltaPhi, deltaEta);

      if (swapped) {
	// swapping the channels
	chan = chan2;
	ieta = ieta2;
	iphi = iphi2;
      }
      
    } // loop on the second digit
  } // loop on the first digit
}

//__________________________________________________
bool Cluster::mergeClusters(Cluster& c, Digit* dig, int deltaPhi, int deltaEta){

  // routine to merge two clusters

  float totDigit = dig->getTOT()*Geo::TOTBIN*1E-3;
  float timeDigit = dig->getTDC()*Geo::TDCBIN;
  int mask;
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
      
  dig->setIsUsedInCluster();

  if (c.getTot() < totDigit) {
    // the main channel contributing to the cluster is the one from the digit --> changing it
    c.setTot(totDigit);
    c.setTime(timeDigit);
    c.setL0L1Latency(/*to be set from the digit*/);
    c.setDeltaBC(/*to be set from the digit*/ );
    c.setMainContributingChannel(dig->getChannel());
    c.addAndShiftContributingChannels(mask);
    return kTRUE;
  }
  else {
    c.addBitInContributingChannels(mask);
  }

  return kFALSE;

}
  
  
  



      

//__________________________________________________
void Clusterer::initStrip()
{
  mPrev = mColumn1 + 1;
  mCurr = mColumn2 + 1;
  std::fill(std::begin(mColumn1), std::end(mColumn1), -1);
  std::fill(std::begin(mColumn2), std::end(mColumn2), -1);

  mDigits.clear();
  mPreClusterHeads.clear();
  mPreClusterIndices.clear();
  Digit* dig = &mStripData.digits[0];
  //mCol = pix->col; // not for TOF
  //mCurr[pix->row] = 0; // not for TOF
  // start the first pre-cluster
  mPreClusterHeads.push_back(0);
  mPreClusterIndices.push_back(0);
  mDigits.emplace_back(-1, dig);
}

//__________________________________________________
void Clusterer::updateStrip(int idig)
{
  Digit* dig = &mStripData.digits[idig];

  // up to here for now...


















  

  if (mCol != pix->col) { // switch the buffers
    Int_t* tmp = mCurr;
    mCurr = mPrev;
    mPrev = tmp;
    if (pix->col > mCol + 1)
      std::fill(mPrev, mPrev + kMaxRow, -1);
    std::fill(mCurr, mCurr + kMaxRow, -1);
    mCol = pix->col;
  }

  Bool_t attached = false;
  UShort_t row = pix->row;
  Int_t neighbours[]{ mCurr[row - 1], mPrev[row], mPrev[row + 1], mPrev[row - 1] };
  for (auto pci : neighbours) {
    if (pci < 0)
      continue;
    auto& ci = mPreClusterIndices[pci];
    if (attached) {
      auto& newci = mPreClusterIndices[mCurr[row]];
      if (ci < newci)
        newci = ci;
      else
        ci = newci;
    } else {
      auto& firstIndex = mPreClusterHeads[ci];
      mPixels.emplace_back(firstIndex, pix);
      firstIndex = mPixels.size() - 1;
      mCurr[row] = pci;
      attached = true;
    }
  }

  if (attached)
    return;

  // start new precluster
  mPreClusterHeads.push_back(mPixels.size());
  mPixels.emplace_back(-1, pix);
  Int_t lastIndex = mPreClusterIndices.size();
  mPreClusterIndices.push_back(lastIndex);
  mCurr[row] = lastIndex;
}

