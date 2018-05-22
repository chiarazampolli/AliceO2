// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//
//  TOF strip class: it will be used to store the hits and digits at TOF that
//  fall in the same strip
//

#ifndef ALICEO2_TOF_STRIP_
#define ALICEO2_TOF_STRIP_

#include <TOF/Digit.h>
#include <TObject.h> // for TObject
#include <exception>
#include <map>
#include <sstream>
#include <vector>
#include "MathUtils/Cartesian3D.h"
#include "SimulationDataFormat/MCCompLabel.h"


namespace o2
{
namespace TOF
{

/// @class Strip
/// @brief Container for similated points connected to a given TOF strip
/// This will be used in order to allow a more efficient clusterization
/// that can happen only between digits that belong to the same strip
///
.
  class Strip{

  public:
    /// Default constructor
  Strip() = default;

  /// Destructor
  ~Strip() = default;

  /// Main constructor
  /// @param stripindex Index of the strip
  /// @param mat Transformation matrix
  Strip(Int_t index);

  /// Copy constructor
  /// @param ref Reference for the copy
  Strip(const Strip& ref);

  /// Empties the point container
  /// @param option unused
  void Clear();

  /// Change the chip index
  /// @param index New chip index
  void SetStripIndex(Int_t index) { mStripIndex = index; }
  void Init(Int_t index, const o2::Transform3D* mat)
  {
    mStripIndex = index;
    mMat = mat;
  }

  /// Get the chip index
  /// @return Index of the chip
  Int_t GetStripIndex() const { return mStripIndex; }
  /// Insert new ITSMFT point into the Chip
  /// @param p Hit to be added
  void InsertHit(const HitType* p);

  /// Get the number of point assigned to the chip
  /// @return Number of points assigned to the chip
  Int_t GetNumberOfHits() const { return mHits.size(); }

  /// Get the strip index from hit
  Int_t GetStripIndex(const HitType* hit);

  /// reset points container
  void ClearHits() { mHits.clear(); }
  o2::TOF::Digit* findDigit(ULong64_t key);

  /// Access Hit assigned to chip at a given index
  /// @param index Index of the point
  /// @return Hit at given index (nullptr if index is out of bounds)
  const HitType* GetHitAt(Int_t index) const;  
  
  Int_t addDigit(UInt_t roframe, UShort_t row, UShort_t col, float charge, Int_t lbl, double timestamp); // returns the MC label 

  void fillOutputContainer(std::vector<Digit>* digits, UInt_t maxFrame);

 protected:
  Int_t mStripIndex = -1;                          ///< Strip ID
  std::vector<const HitType*> mHits;                  ///< Hits connnected to the given chip
  std::map<ULong64_t, o2::TOF::Digit> mDigits; ///< Map of fired pixels, possibly in multiple frames

  ClassDefNV(Strip, 1);
};

inline o2::TOF::Digit* Strip::findDigit(ULong64_t key)
{
  // finds the digit corresponding to global key
  auto digitentry = mDigits.find(key);
  return digitentry != mDigits.end() ? &(digitentry->second) : nullptr;
}


#endif /* defined(ALICEO2_TOF_STRIP_) */

  
