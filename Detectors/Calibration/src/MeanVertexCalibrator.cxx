// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DetectorsCalibration/MeanVertexCalibrator.h"
#include "Framework/Logger.h"
#include "MathUtils/MathBase.h"
#include "CommonUtils/MemFileHelper.h"
#include "CCDB/CcdbApi.h"
#include "DetectorsCalibration/Utils.h"

namespace o2
{
namespace calibration
{

using Slot = o2::calibration::TimeSlot<o2::calibration::MeanVertexData>;
using o2::math_utils::math_base::fitGaus;
using clbUtils = o2::calibration::Utils;
  using MeanVertexObject = o2::dataformats::MeanVertexObject;

void MeanVertexCalibrator::initOutput()
{
  // Here we initialize the vector of our output objects
  mInfoVector.clear();
  mMeanVertexVector.clear();
  return;
}

//_____________________________________________
void MeanVertexCalibrator::finalizeSlot(Slot& slot)
{
  // Extract results for the single slot
  o2::calibration::MeanVertexData* c = slot.getContainer();
  LOG(INFO) << "Finalize slot " << slot.getTFStart() << " <= TF <= " << slot.getTFEnd() << " with "
            << c->getEntries() << " entries";
  MeanVertexObject mvo; 
  if (mUseFit) {
    // x coordinate
    std::vector<float> fitValues;
    float* array = &c->histoX[0];
    int fitres = fitGaus(c->nbinsX, array, -(c->rangeX), c->rangeX, fitValues);
    if (fitres >= 0) {
      LOG(INFO) << "X: Fit result " << fitres << " Mean = " << fitValues[1] << " Sigma = " << fitValues[2];
    } else {
      LOG(ERROR) << "X: Fit failed with result = " << fitres;
    }
    mvo.setX(fitValues[1]);
    mvo.setSigmaX(fitValues[2]);

    // y coordinate
    array = &c->histoY[0];
    fitres = fitGaus(c->nbinsY, array, -(c->rangeY), c->rangeY, fitValues);
    if (fitres >= 0) {
      LOG(INFO) << "Y: Fit result " << fitres << " Mean = " << fitValues[1] << " Sigma = " << fitValues[2];
    } else {
      LOG(ERROR) << "Y: Fit failed with result = " << fitres;
    }
    mvo.setY(fitValues[1]);
    mvo.setSigmaY(fitValues[2]);

    // z coordinate
    array = &c->histoZ[0];
    fitres = fitGaus(c->nbinsZ, array, -(c->rangeZ), c->rangeZ, fitValues);
    if (fitres >= 0) {
      LOG(INFO) << "Z: Fit result " << fitres << " Mean = " << fitValues[1] << " Sigma = " << fitValues[2];
    } else {
      LOG(ERROR) << "Z: Fit failed with result = " << fitres;
    }
    mvo.setZ(fitValues[1]);
    mvo.setSigmaZ(fitValues[2]);
    mTmpMVVector.emplace(mTmpMVVector.begin() + slot.getTFStart(), std::move(mvo));
  }
  else {
    // for now I do nothing
    //LOG(FATAL) << "For now, only fit is allowed";
    mTmpMVData.emplace(mTmpMVData.begin() + slot.getTFStart(), std::move(c));
  }

  // TODO: the timestamp is now given with the TF index, but it will have
  // to become an absolute time. This is true both for the lhc phase object itself
  // and the CCDB entry
  std::map<std::string, std::string> md;
  auto clName = o2::utils::MemFileHelper::getClassName(mvo);
  auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
  mInfoVector.emplace_back("GRP/MeanVertex", clName, flName, md, slot.getTFStart(), 99999999999999);
  mMeanVertexVector.emplace_back(mvo);

  slot.print();
}

//_____________________________________________
Slot& MeanVertexCalibrator::emplaceNewSlot(bool front, TFType tstart, TFType tend)
{
  auto& cont = getSlots();
  auto& slot = front ? cont.emplace_front(tstart, tend) : cont.emplace_back(tstart, tend);
  slot.setContainer(std::make_unique<MeanVertexData>(mUseFit, mNBinsX, mRangeX, mNBinsY, mRangeY, mNBinsZ, mRangeZ));
  return slot;
}

} // end namespace calibration
} // end namespace o2
