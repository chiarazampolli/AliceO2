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
    for (int idig = 1; idig < mStripData.digits.size(); idig++)
      updateStrip(idig);
    finishStrip(clusters);
  }
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

