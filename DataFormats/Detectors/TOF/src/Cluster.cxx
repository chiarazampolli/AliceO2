// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Cluster.cxx
/// \brief Implementation of the TOF cluster

#include "DataFormatsTOF/Cluster.h"
#include "FairLogger.h"
#include "TOFBase/Geo.h"

#include <TMath.h>
#include <TString.h>

#include <cstdlib>

using namespace o2::TOF;

ClassImp(o2::TOF::Cluster);

Cluster::Cluster(std::int16_t sensid, float x, float y, float z, float sy2, float sz2, float syz, float timeRaw, float time,float tot, int L0L1Latency, int deltaBC) : o2::BaseCluster<float>(sensid, x, y, z, sy2, sz2, syz), mTimeRaw(timeRaw), mTime(time), mTot(tot), mL0L1Latency(L0L1Latency), mDeltaBC(deltaBC) {

  // caching R and phi
  mR = TMath::Sqrt(x*x + y*y);
  mPhi = TMath::ATan2(y, x);
  mContributingChannels = 0;
  
}

//______________________________________________________________________
int Cluster::getNumOfContributingChannels() const {
  //
  // returning how many hits contribute to this cluster
  //
  int nContributingChannels = 0;
  if (mContributingChannels == 0) {
    LOG(ERROR) << "The current cluster has no hit contributing to it!" << FairLogger::endl;
  }
  else {
    nContributingChannels++;
    if ((mContributingChannels & 0x40000) == 0x40000)  nContributingChannels++; // alsoUP
    if ((mContributingChannels & 0x80000) == 0x80000)  nContributingChannels++; // alsoDOWN
    if ((mContributingChannels & 0x100000) == 0x100000)  nContributingChannels++; // alsoRIGHT
    if ((mContributingChannels & 0x200000) == 0x200000)  nContributingChannels++; // alsoLEFT
  }
  return nContributingChannels;
}  

//______________________________________________________________________
std::ostream& operator<<(std::ostream& os, const Cluster& c) {  
  os << (o2::BaseCluster<float>&)c;
  os << " TOF cluster: raw time = " << std::scientific << c.getTimeRaw() << ", time = "  << std::scientific << c.getTime() << ", Tot = " << std::scientific << c.getTot() << ", L0L1Latency = " << c.getL0L1Latency() << ", deltaBC = " << c.getDeltaBC() << ", R = " << c.getR() << ", mPhi = " << c.getPhi() << ", ContributingChannels = " << c.getNumOfContributingChannels() << "\n";
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
void Cluster::buildCluster() {

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

  setMainContributingChannel(mContributingDigit[0]->getChannel());
  setTime(mContributingDigit[0]->getTDC()*Geo::TDCBIN); // time in ps (for now we assume it calibrated)
  setTot(mContributingDigit[0]->getTOT()*Geo::TOTBIN*1E-3); // TOT in ns (for now we assume it calibrated)
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
    addBitInContributingChannels(mask);
  }

  return;
  
}
    
    
