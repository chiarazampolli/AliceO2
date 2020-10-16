// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef MEAN_VERTEX_DATA_H_
#define MEAN_VERTEX_DATA_H_

#include "DetectorsCalibration/TimeSlot.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include <array>
#include <gsl/span>

namespace o2
{
namespace calibration
{

struct MeanVertexData {

  using PVertex = o2::dataformats::PrimaryVertex;
  float rangeX = 10.f;
  float rangeY = 10.f;
  float rangeZ = 10.f;
  int nbinsX = 1000;
  int nbinsY = 1000;
  int nbinsZ = 1000;
  float v2BinX = nbinsX / (2 * rangeX);
  float v2BinY = nbinsY / (2 * rangeY);
  float v2BinZ = nbinsZ / (2 * rangeZ);
  int entries = 0;
  std::vector<float> histoX{0};
  std::vector<float> histoY{0};
  std::vector<float> histoZ{0};
  bool useFit = false;
  
  MeanVertexData();

MeanVertexData(bool buseFit, int nbX, float rX, int nbY, float rY, int nbZ, float rZ) :
  useFit(buseFit), nbinsX(nbX), rangeX(rX), v2BinX(0), nbinsY(nbY), rangeY(rY), v2BinY(0),
    nbinsZ(nbZ), rangeZ(rZ), v2BinZ(0) 
  {
    if (useFit) {
      if (rX <= 0. || nbX < 1 || rY <= 0. || nbY < 1 || rZ <= 0. || nbZ < 1) {
	throw std::runtime_error("Wrong initialization of the histogram");
      }
      v2BinX = nbinsX / (2 * rangeX);
      v2BinY = nbinsY / (2 * rangeY);
      v2BinZ = nbinsZ / (2 * rangeZ);
      histoX.resize(nbinsX, 0.);
      histoY.resize(nbinsY, 0.);
      histoZ.resize(nbinsZ, 0.);
    }
  }

  size_t getEntries() const { return entries; }
  void print() const;
  void fill(const gsl::span<const PVertex> data);
  void merge(const MeanVertexData* prev);

  ClassDefNV(MeanVertexData, 1);
};

} // end namespace calibration
} // end namespace o2

#endif /* MEAN_VERTEX_DATA_H_ */
