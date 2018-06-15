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
void Cluster::addAndShiftContributingChannels(int mask) {

  // changing the existing contributing channels wrt the new main one
  // and setting the new main one

  if (mask > kRight)
    mask = mask << 4; // move to the left
  else
    mask = mask >> 4; // move to the right
  addBitInContributingChannels(mask);

  // to be done to move the other bits

  return;
}

