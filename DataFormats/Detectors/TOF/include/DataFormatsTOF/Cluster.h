// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Cluster.h
/// \brief Definition of the TOF cluster

#ifndef ALICEO2_TOF_CLUSTER_H
#define ALICEO2_TOF_CLUSTER_H

#include "ReconstructionDataFormats/BaseCluster.h"
#include <boost/serialization/base_object.hpp> // for base_object


namespace o2
{
namespace TOF
{
/// \class Cluster
/// \brief Cluster class for TOF
///

class Cluster : public o2::BaseCluster<float>
{
 public:
  Cluster() = default;

  Cluster(std::int16_t sensid, float x, float y, float z, float sy2, float sz2, float syz, float timeRaw, float time, float tot, int L0L1latency, int deltaBC);
  
  ~Cluster() = default;

  float  getTimeRaw()    const {return mTimeRaw;}    // Cluster ToF getter
  void   setTimeRaw(float timeRaw)   {mTimeRaw = timeRaw;}       // Cluster ToF setter
  float  getTime()    const {return mTime;}    // Cluster ToF getter
  void   setTime(float time)   {mTime = time;}       // Cluster ToF setter
  float  getTot()    const {return mTot;}    // Cluster Charge getter
  void   setTot(int tot)       {mTot = tot;}       // Cluster ToT setter  
  int    getL0L1Latency() const    {return mL0L1Latency;}; // L0L1 latency
  void   setL0L1Latency(int value) {mL0L1Latency = value;}; // L0-L1 latency
  int    getDeltaBC() const    {return mDeltaBC;}; // deltaBC
  void   setDeltaBC(int value) {mDeltaBC = value;}; // deltaBC
  //float  getZ()   const   {return mZ;}   // Cluster Z - already in the definition of the cluster  
  float  getR()   const   {return mR;}   // Cluster Radius  
  float  getPhi() const   {return mPhi;} // Cluster Phi

  int    getContributingChannels() const {return mContributingChannels;}
  int    getNumOfContributingChannels() const; // returns the number of hits associated to the cluster, i.e. the number of hits that built the cluster; it is the equivalente of teh old AliESDTOFCluster::GetNTOFhits() 
  int    getMainContributingChannel() const {return mContributingChannels & 0x3FFFF;}
  bool   isUpContributing() const {return (mContributingChannels & 0x40000) == 0x40000 ? 1 : 0;}
  bool   isDownContributing() const {return (mContributingChannels & 0x80000) == 0x80000 ? 1 : 0;}
  bool   isLeftContributing() const {return (mContributingChannels & 0x100000) == 0x100000 ? 1 : 0;}
  bool   isRightContributing() const {return (mContributingChannels & 0x200000) == 0x200000 ? 1 : 0;}
  
 private:
  friend class boost::serialization::access;

  float  mTimeRaw;     // raw TOF time // CZ: in AliRoot it is a double
  float  mTime;        // calibrated TOF time // CZ: in AliRoot it is a double
  float  mTot;         // Time-Over-threshold // CZ: in AliRoot it is a double
  int    mL0L1Latency; // L0L1 latency // CZ: is it different per cluster? Checking one ESD file, it seems that it is always the same (see: /alice/data/2017/LHC17n/000280235/pass1/17000280235019.100/AliESDs.root)
  int    mDeltaBC;     // DeltaBC --> can it be a char or short? // CZ: is it different per cluster? Checking one ESD file, it seems that it can vary (see: /alice/data/2017/LHC17n/000280235/pass1/17000280235019.100/AliESDs.root)
  //float  mZ;           //! z-coordinate // CZ: to be verified if it is the same in the BaseCluster class
  float  mR;           //! radius
  float  mPhi;         //! phi coordinate 
  int    mContributingChannels;       // index of the channels that contributed to the cluster; to be read like this:
                                      //channel & 0x3FFFF -> first 18 bits to strore the main channel
                                      //channel & bit19 (0x40000) -> alsoUP
                                      //channel & bit20 (0x80000)-> alsoDOWN
                                      //channel & bit21 (0x100000)-> alsoRIGHT
                                      //channel & bit22 (0x200000)-> alsoLEFT
  
  ClassDefNV(Cluster, 1);
};

 std::ostream& operator<<(std::ostream& os, const Cluster& c);
} // namespace TOF
} // namespace o2
#endif
