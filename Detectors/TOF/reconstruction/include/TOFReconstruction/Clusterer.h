// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Clusterer.h
/// \brief Definition of the TOF cluster finder
#ifndef ALICEO2_TOF_CLUSTERER_H
#define ALICEO2_TOF_CLUSTERER_H

#include <utility>
#include <vector>
#include "DataFormatsTOF/Cluster.h"
#include "TOFBase/Geo.h"
#include "TOFReconstruction/DataReader.h"

namespace o2
{
class MCCompLabel;
namespace dataformats
{
template <typename T>
class MCTruthContainer;
}

namespace tof
{
class Clusterer
{
  //using Cluster = o2::TOF::Cluster;
  using Label = o2::MCCompLabel;

 public:
  Clusterer();
  ~Clusterer() = default;
  
  Clusterer(const Clusterer&) = delete;
  Clusterer& operator=(const Clusterer&) = delete;

  void process(DataReader& r, std::vector<Cluster>& clusters);

  void setMCTruthContainer(o2::dataformats::MCTruthContainer<o2::MCCompLabel>* truth) { mClsLabels = truth; }

 private:

  void initStrip();
  void updateStrip(int ip);
  void finishStrip(std::vector<Cluster>& clusters);
  void fetchMCLabels(const Digit* dig, std::array<Label, Cluster::maxLabels>& labels, int& nfilled) const;

  StripData mStripData; ///< single strip data provided by the reader

  Int_t mColumn1[kMaxRow + 2]; // not sure we will need it, for now I just took it from ITSMFT
  Int_t mColumn2[kMaxRow + 2]; // not sure we will need it, for now I just took it from ITSMFT
  Int_t *mCurr, *mPrev; // not sure we will need it, for now I just took it from ITSMFT

  using NextIndex = Int_t;
  std::vector<std::pair<NextIndex, const Digit*>> mDigits;

  using FirstIndex = Int_t;
  std::vector<FirstIndex> mPreClusterHeads;

  std::vector<Int_t> mPreClusterIndices;

  UShort_t mCol = 0xffff; ///< Column being processed

  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mClsLabels = nullptr; // Cluster MC labels
};

} // namespace tof
} // namespace o2
#endif /* ALICEO2_TOF_CLUSTERER_H */
