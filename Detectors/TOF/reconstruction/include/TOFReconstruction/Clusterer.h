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

namespace o2
{
class MCCompLabel;
namespace dataformats
{
template <typename T>
class MCTruthContainer;
}

namespace TOF
{
class Clusterer
{
  using Cluster = o2::TOF::Cluster;
  using Label = o2::MCCompLabel;

 public:
  Clusterer();
  ~Clusterer() = default;
  
  Clusterer(const Clusterer&) = delete;
  Clusterer& operator=(const Clusterer&) = delete;

  void process(PixelReader& r, std::vector<Cluster>& clusters);

  // provide the common ITSMFT::GeometryTGeo to access matrices
  void setGeometry(const o2::TOF::Geo* gm) { mGeometry = gm; }
  void setMCTruthContainer(o2::dataformats::MCTruthContainer<o2::MCCompLabel>* truth) { mClsLabels = truth; }

 private:

  ...

   const o2::TOF::Geo* mGeometry = nullptr;                      ///< TOF upgrade geometry
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* mClsLabels = nullptr; // Cluster MC labels
};

} // namespace TOF
} // namespace o2
#endif /* ALICEO2_TOF_CLUSTERER_H */
