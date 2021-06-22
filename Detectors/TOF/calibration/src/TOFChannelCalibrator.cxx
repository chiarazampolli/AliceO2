// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TOFCalibration/TOFChannelCalibrator.h"
#include "Framework/Logger.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <TStopwatch.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace o2
{
namespace tof
{

using Slot = o2::calibration::TimeSlot<o2::tof::TOFChannelData>;
using TimeSlewing = o2::dataformats::CalibTimeSlewingParamTOF;
using clbUtils = o2::calibration::Utils;
using boost::histogram::indexed;
using namespace o2::tof;
//using boost::histogram::algorithm; // not sure why it does not work...

//_____________________________________________
void TOFChannelData::fill(const gsl::span<const o2::dataformats::CalibInfoTOF> data)
{
  // fill container
  for (int i = data.size(); i--;) {
    auto ch = data[i].getTOFChIndex();
    int sector = ch / Geo::NPADSXSECTOR;
    int chInSect = ch % Geo::NPADSXSECTOR;
    auto dt = data[i].getDeltaTimePi();

    auto tot = data[i].getTot();
    // TO BE DISCUSSED: could it be that the LHCphase is too old? If we ar ein sync mode, it could be that it is not yet created for the current run, so the one from the previous run (which could be very old) is used. But maybe it does not matter much, since soon enough a calibrated LHC phase should be produced
    auto corr = mCalibTOFapi->getTimeCalibration(ch, tot); // we take into account LHCphase, offsets and time slewing
    LOG(DEBUG) << "inserting in channel " << ch << ": dt = " << dt << ", tot = " << tot << ", corr = " << corr << ", corrected dt = " << dt - corr;

    dt -= corr;

    mHisto[sector](dt, chInSect); // we pass the calibrated time
    mEntries[ch] += 1;
  }
}

//_____________________________________________
void TOFChannelData::fill(const gsl::span<const o2::tof::CalibInfoCluster> data)
{
  // fill container
  for (int i = data.size(); i--;) {
    auto ch = data[i].getCH();
    auto dch = data[i].getDCH(); // this is a char! if you print it, you need to cast it to int
    auto dt = data[i].getDT();
    auto tot1 = data[i].getTOT1();
    auto tot2 = data[i].getTOT2();

    // we order them so that the channel number of the first cluster is smaller than
    // the one of the second cluster
    if (dch < 0) {
      ch += dch;
      dt = -dt;
      dch = -dch;
      float inv = tot1;
      tot1 = tot2;
      tot2 = inv;
    }

    int sector = ch / Geo::NPADSXSECTOR;
    int absoluteStrip = ch / Geo::NPADS;
    int stripInSect = absoluteStrip % Geo::NSTRIPXSECTOR;
    int shift = 0;
    if (dch == 1) {
      shift = 0; // 2nd channel is on the right
    } else if (dch == 48) {
      shift = 1; // 2nd channel is at the top
    } else {
      continue;
    }
    int chOnStrip = ch % 96;
    int comb = 96 * shift + chOnStrip; // index of the current pair of clusters on the strip; there are in total 96 + 48

    auto corr1 = mCalibTOFapi->getTimeCalibration(ch, tot1);       // we take into account LHCphase, offsets and time slewing
    auto corr2 = mCalibTOFapi->getTimeCalibration(ch + dch, tot2); // we take into account LHCphase, offsets and time slewing
    LOG(DEBUG) << "inserting in channel " << ch << ", " << ch + dch << ": dt = " << dt << ", tot1 = " << tot1 << ", tot2 = " << tot2 << ", corr1 = " << corr1 << ", corr2 = " << corr2 << ", corrected dt = " << dt - corr1 + corr2;

    dt -= corr1 - corr2;

    int combInSect = comb + stripInSect * NCOMBINSTRIP;

    LOG(DEBUG) << "ch = " << ch << ", sector = " << sector << ", absoluteStrip = " << absoluteStrip << ", stripInSect = " << stripInSect << ", shift = " << shift << ", dch = " << (int)dch << ", chOnStrip = " << chOnStrip << ", comb = " << comb;

    mHisto[sector](dt, combInSect); // we pass the difference of the *calibrated* times
    mEntries[comb + NCOMBINSTRIP * absoluteStrip] += 1;
  }
}

//_____________________________________________
void TOFChannelData::merge(const TOFChannelData* prev)
{
  // merge data of 2 slots
  for (int isect = 0; isect < Geo::NSECTORS; isect++) {
    mHisto[isect] += prev->getHisto(isect);
  }
  for (auto iel = 0; iel < mEntries.size(); iel++) {
    mEntries[iel] += prev->mEntries[iel];
  }
}

//_____________________________________________
bool TOFChannelData::hasEnoughData(int minEntries) const
{
  // true if all channels can be fitted --> have enough statistics

  // We consider that we have enough entries if the mean of the number of entries in the channels
  // with at least one entry is greater than the cut at "minEntries"
  // Channels/pairs with zero entries are assumed to be off --> we do not consider them

  //printEntries();
  int nValid = 0;
  float mean = 0;
  int smallestElementIndex = -1;
  int smallestEntries = 1e5;
  int largestElementIndex = -1;
  int largestEntries = 0;
  for (auto i = 0; i < mEntries.size(); ++i) {
    if (mEntries[i] != 0) { // skipping channels/pairs if they have zero entries (most likely they are simply off)
      mean += mEntries[i];
      ++nValid;
      if (mEntries[i] < smallestEntries) {
        smallestEntries = mEntries[i];
        smallestElementIndex = i;
      }
      if (mEntries[i] > largestEntries) {
        largestEntries = mEntries[i];
        largestElementIndex = i;
      }
    }
  }
  if (nValid == 0) {
    LOG(INFO) << "hasEnough = false: all channels/pairs are empty";
    return false;
  }

  LOG(DEBUG) << "mean = " << mean << ", nvalid = " << nValid;
  mean /= nValid;

  LOG(DEBUG) << "minElement is at position " << smallestElementIndex << " and has " << smallestEntries << " entries";
  LOG(DEBUG) << "maxElement is at position " << largestElementIndex << " and has " << largestEntries << " entries";
  float threshold = minEntries + 5 * std::sqrt(minEntries);
  bool enough = mean < threshold ? false : true;
  if (enough) {
    LOG(INFO) << "hasEnough: " << (int)enough << " ("
              << nValid << " valid channels found (should be " << mEntries.size() << ") with mean = "
              << mean << " with cut at = " << threshold << ") ";
  }
  return enough;
}

//_____________________________________________
void TOFChannelData::print() const
{
  LOG(INFO) << "Printing histograms:";
  std::ostringstream os;
  for (int isect = 0; isect < Geo::NSECTORS; isect++) {
    LOG(INFO) << "Sector: " << isect;
    os << mHisto[isect];
    auto nentriesInSec = boost::histogram::algorithm::sum(mHisto[isect]);
    LOG(INFO) << "Number of entries in histogram: " << boost::histogram::algorithm::sum(mHisto[isect]);
    int cnt = 0;
    if (nentriesInSec != 0) {
      for (auto&& x : indexed(mHisto[isect])) {
        if (x.get() > 0) {
          const auto i = x.index(0); // current index along first axis --> t-texp
          const auto j = x.index(1); // current index along second axis --> channel
          const auto b0 = x.bin(0);  // current bin interval along first axis --> t-texp
          const auto b1 = x.bin(1);  // current bin interval along second axis --> channel
          LOG(INFO) << "bin " << cnt << ": channel = " << j << " in [" << b1.lower() << ", " << b1.upper()
                    << "], t-texp in [" << b0.lower() << ", " << b0.upper() << "], has entries = " << x.get();
        }
        cnt++;
      }
    }
    LOG(INFO) << cnt << " bins inspected";
  }
}

//_____________________________________________
void TOFChannelData::print(int isect) const
{
  LOG(INFO) << "*** Printing histogram " << isect;
  std::ostringstream os;
  int cnt = 0;
  os << mHisto[isect];
  LOG(INFO) << "Number of entries in histogram: " << boost::histogram::algorithm::sum(mHisto[isect]);
  for (auto&& x : indexed(mHisto[isect])) { // does not work also when I use indexed(*(mHisto[sector]))
    cnt++;
    LOG(DEBUG) << " c " << cnt << " i " << x.index(0) << " j " << x.index(1) << " b0 " << x.bin(0) << " b1 " << x.bin(1) << " val= " << *x << "|" << x.get();
    if (x.get() > 0) {
      LOG(INFO) << "x = " << x.get() << " c " << cnt;
    }
  }
  LOG(INFO) << cnt << " bins inspected";
}

//_____________________________________________
void TOFChannelData::printEntries() const
{
  // to print number of entries per channel
  for (int i = 0; i < mEntries.size(); ++i) {
    if (mEntries.size() > tof::Geo::NCHANNELS) {
      LOG(INFO) << "pair of channels " << i << " has " << mEntries[i] << " entries";
    } else {
      LOG(INFO) << "channel " << i << " has " << mEntries[i] << " entries";
    }
  }
}

//_____________________________________________
int TOFChannelData::findBin(float v) const
{
  // find the bin along the x-axis (with t-texp) where the value "v" is; this does not depend on the channel
  // (axis 1), nor on the sector, so we use sector0

  if (v == mRange) {
    v -= 1.e-1;
  }

  LOG(DEBUG) << "In FindBin, v = : " << v;
  LOG(DEBUG) << "bin0 limits: lower = " << mHisto[0].axis(0).bin(0).lower() << ", upper = " << mHisto[0].axis(0).bin(0).upper();
  LOG(DEBUG) << "bin1000 limits: lower = " << mHisto[0].axis(0).bin(mNBins - 1).lower() << ", upper = " << mHisto[0].axis(0).bin(mNBins - 1).upper();
  LOG(DEBUG) << "v = " << v << " is in bin " << mHisto[0].axis(0).index(v);

  return mHisto[0].axis(0).index(v);
}

//_____________________________________________
float TOFChannelData::integral(int chmin, int chmax, float binmin, float binmax) const
{
  // calculates the integral in [chmin, chmax] and in [binmin, binmax]

  if (binmin < -mRange || binmax > mRange || chmin < 0 || chmax >= Geo::NSECTORS * mNElsPerSector) {
    throw std::runtime_error("Check your bins, we cannot calculate the integrals in under/overflows bins");
  }
  if (binmax < binmin || chmax < chmin) {
    throw std::runtime_error("Check your bin limits!");
  }

  int sector = chmin / mNElsPerSector;
  if (sector != chmax / mNElsPerSector) {
    throw std::runtime_error("We cannot integrate over channels that belong to different sectors");
  }

  int chinsectormin = chmin % mNElsPerSector;
  int chinsectormax = chmax % mNElsPerSector;

  float res2 = 0;
  //TStopwatch t3;
  int ind = -1;
  int binxmin = findBin(binmin);
  int binxmax = findBin(binmax);
  LOG(DEBUG) << "binxmin = " << binxmin << ", binxmax = " << binxmax;
  //t3.Start();
  for (unsigned j = chinsectormin; j <= chinsectormax; ++j) {
    for (unsigned i = binxmin; i <= binxmax; ++i) {
      const auto& v = mHisto[sector].at(i, j);
      res2 += v;
    }
  }
  //t3.Stop();
  LOG(DEBUG) << "Time for integral looping over axis (result = " << res2 << "):";
  //t3.Print();

  return res2;

}

//_____________________________________________
float TOFChannelData::integral(int ch, float binmin, float binmax) const
{
  // calculates the integral along one fixed channel and in [binmin, binmax]

  return integral(ch, ch, binmin, binmax);
}

//_____________________________________________
float TOFChannelData::integral(int chmin, int chmax, int binxmin, int binxmax) const
{
  // calculates the integral in [chmin, chmax] and in [binmin, binmax]

  if (binxmin < 0 || binxmax > mNBins || chmin < 0 || chmax >= Geo::NSECTORS * mNElsPerSector) {
    throw std::runtime_error("Check your bins, we cannot calculate the integrals in under/overflows bins");
  }
  if (binxmax < binxmin || chmax < chmin) {
    throw std::runtime_error("Check your bin limits!");
  }

  int sector = chmin / mNElsPerSector;
  if (sector != chmax / mNElsPerSector) {
    throw std::runtime_error("We cannot integrate over channels that belong to different sectors");
  }

  int chinsectormin = chmin % mNElsPerSector;
  int chinsectormax = chmax % mNElsPerSector;

  float res2 = 0;
  //TStopwatch t3;
  //t3.Start();
  for (unsigned j = chinsectormin; j <= chinsectormax; ++j) {
    for (unsigned i = binxmin; i <= binxmax; ++i) {
      const auto& v = mHisto[sector].at(i, j);
      res2 += v;
    }
  }
  //t3.Stop();
  LOG(DEBUG) << "Time for integral looping over axis (result = " << res2 << "):";
  //t3.Print();
  return res2;

}

//_____________________________________________
float TOFChannelData::integral(int ch, int binxmin, int binxmax) const
{
  // calculates the integral along one fixed channel and in [binmin, binmax]

  return integral(ch, ch, binxmin, binxmax);
}

//_____________________________________________
float TOFChannelData::integral(int ch) const
{
  // calculates the integral along one fixed channel and in the full x-range

  return integral(ch, ch, 0, mNBins - 1);
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

template <typename T>
void TOFChannelCalibrator<T>::finalizeSlotWithCosmics(Slot& slot)
{    // Extract results for the single slot
    o2::tof::TOFChannelData* c = slot.getContainer();
    LOG(INFO) << "Finalize slot for calibration with cosmics " << slot.getTFStart() << " <= TF <= " << slot.getTFEnd();
    LOG(INFO) << "Finalize slot for calibration with cosmics 1) " << slot.getTFStart() << " <= TF <= " << slot.getTFEnd();

    if (!mFuncDeltaOffset) {
      mFuncDeltaOffset = new TF1("fChanOffset", FuncDeltaOffset, 0, NCOMBINSTRIP, Geo::NPADS);
      mFuncDeltaOffset->FixParameter(24, 0); // fix one pad to fixed offset (as reference)
      for (int ichLocal = 0; ichLocal < Geo::NPADS; ichLocal++) {
        if (ichLocal == 24) {
          continue; //already fixed
          mFuncDeltaOffset->SetParLimits(ichLocal, -mRange, mRange);
        }
      }
    }
    LOG(INFO) << "Finalize slot for calibration with cosmics 2) " << slot.getTFStart() << " <= TF <= " << slot.getTFEnd();

    // for the CCDB entry
    std::map<std::string, std::string> md;
    TimeSlewing& ts = mCalibTOFapi->getSlewParamObj(); // we take the current CCDB object, since we want to simply update the offset

    float xp[NCOMBINSTRIP], exp[NCOMBINSTRIP], deltat[NCOMBINSTRIP], edeltat[NCOMBINSTRIP], fracUnderPeak[Geo::NPADS];
    LOG(INFO) << "Number of threads that will be used 0) = " << mNThreads;

#ifdef WITH_OPENMP
    LOG(INFO) << "omp_get_max_threads() = " << omp_get_max_threads();
    if (mNThreads < 1) mNThreads = omp_get_max_threads();
    LOG(INFO) << "Number of threads that will be used = " << mNThreads;
#pragma omp parallel for schedule(dynamic) num_threads(mNThreads)
#endif
    for (int sector = 0; sector < Geo::NSECTORS; sector++) {
      LOG(INFO) << "Processing sector " << sector;
      int offsetsector = sector * Geo::NSTRIPXSECTOR * Geo::NPADS;
      for (int istrip = 0; istrip < Geo::NSTRIPXSECTOR; istrip++) {
        int offsetstrip = istrip * Geo::NPADS + offsetsector;
        int goodpoints = 0;

        memset(&fracUnderPeak[0], 0, sizeof(fracUnderPeak));

        for (int ipair = 0; ipair < NCOMBINSTRIP; ipair++) {
          int chinsector = ipair + istrip * NCOMBINSTRIP;
          int ich = chinsector + sector * Geo::NSTRIPXSECTOR * NCOMBINSTRIP;
          auto entriesInPair = c->integral(ich);
          if (entriesInPair < mMinEntries && entriesInPair != 0) {
            LOG(DEBUG) << "pair " << ich << " will not be calibrated since it has only " << entriesInPair << " entries (min = " << mMinEntries << ")";
            LOG(INFO) << "pair " << ich << " will not be calibrated since it has only " << entriesInPair << " entries (min = " << mMinEntries << ")";
            continue;
          }
          // make the slice of the 2D histogram so that we have the 1D of the current channel
          std::vector<float> fitValues;
          std::vector<float> histoValues;
          std::vector<int> entriesPerChannel = c->getEntriesPerChannel();
          if (entriesPerChannel.at(ich) == 0) {
            continue; // skip always since a channel with 0 entries is normal, it will be flagged as problematic
            if (mTest) {
              LOG(INFO) << "Skipping channel " << ich << " because it has zero entries, but it should not be"; // should become error!
              continue;
            } else {
              throw std::runtime_error("We found one channel with no entries, we cannot calibrate!");
            }
          }

          // more efficient way
          auto histo = c->getHisto(sector);
          for (unsigned j = chinsector; j <= chinsector; ++j) {
            for (unsigned i = 0; i < c->getNbins(); ++i) {
              const auto& v = histo.at(i, j);
              LOG(DEBUG) << "channel = " << ich << ", in sector = " << sector << " (where it is channel = " << chinsector << ") bin = " << i << " value = " << v;
              histoValues.push_back(v);
            }
          }

          double fitres = fitGaus(c->getNbins(), histoValues.data(), -(c->getRange()), c->getRange(), fitValues);

          if (fitres >= 0) {
            LOG(DEBUG) << "Pair " << ich << " :: Fit result " << fitres << " Mean = " << fitValues[1] << " Sigma = " << fitValues[2];
          } else {
            LOG(INFO) << "Pair " << ich << " :: Fit failed with result = " << fitres;
            continue;
          }

          if (fitValues[2] < 0) {
            fitValues[2] = -fitValues[2];
          }

          float intmin = fitValues[1] - 5 * fitValues[2]; // mean - 5*sigma
          float intmax = fitValues[1] + 5 * fitValues[2]; // mean + 5*sigma

          if (intmin < -mRange) {
            intmin = -mRange;
          }
          if (intmax < -mRange) {
            intmax = -mRange;
          }
          if (intmin > mRange) {
            intmin = mRange;
          }
          if (intmax > mRange) {
            intmax = mRange;
          }

          xp[goodpoints] = ipair + 0.5;      // pair index
          exp[goodpoints] = 0.0;             // error on pair index (dummy since it is on the pair index)
          deltat[goodpoints] = fitValues[1]; // delta between offsets from channels in pair (from the fit) - in ps
          edeltat[goodpoints] = 20;          // TODO: for now put by default to 20 ps since it was seen to be reasonable; but it should come from the fit: who gives us the error from the fit ??????
          goodpoints++;
          int ch1 = ipair % 96;
          int ch2 = ipair / 96 ? ch1 + 48 : ch1 + 1;
          float fractionUnderPeak = entriesInPair > 0 ? c->integral(ich, intmin, intmax) / entriesInPair : 0;
          // we keep as fractionUnderPeak of the channel the largest one that is found in the 3 possible pairs with that channel (for both channels ch1 and ch2 in the pair)
          if (fracUnderPeak[ch1] < fractionUnderPeak) {
            fracUnderPeak[ch1] = fractionUnderPeak;
          }
          if (fracUnderPeak[ch2] < fractionUnderPeak) {
            fracUnderPeak[ch2] = fractionUnderPeak;
          }
        } // end loop pairs

        // fit strip offset
        if (goodpoints == 0) {
          //LOG(INFO) << "We did not find any good point for strip " << istrip << " in sector " << sector;
          continue;
        }
        LOG(DEBUG) << "We found " << goodpoints << " good points for strip " << istrip << " in sector " << sector << " --> we can fit the TGraph";
        TGraphErrors g(goodpoints, xp, deltat, exp, edeltat);
        g.Fit(mFuncDeltaOffset, "Q0");

        //update calibrations
        for (int ichLocal = 0; ichLocal < Geo::NPADS; ichLocal++) {
          int ich = ichLocal + offsetstrip;
          ts.updateOffsetInfo(ich, mFuncDeltaOffset->GetParameter(ichLocal));
          ts.setFractionUnderPeak(ich / Geo::NPADSXSECTOR, ich % Geo::NPADSXSECTOR, fracUnderPeak[ichLocal]);
          ts.setSigmaPeak(ich / Geo::NPADSXSECTOR, ich % Geo::NPADSXSECTOR, abs(mFuncDeltaOffset->GetParError(ichLocal)));
        }

      } // end loop strips
    }   // end loop sectors

    auto clName = o2::utils::MemFileHelper::getClassName(ts);
    auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
    mInfoVector.emplace_back("TOF/ChannelCalib", clName, flName, md, slot.getTFStart(), 99999999999999);
    mTimeSlewingVector.emplace_back(ts);
}

//_____________________________________________

template <typename T>
void TOFChannelCalibrator<T>::finalizeSlotWithTracks(Slot& slot)
  {
    // Extract results for the single slot
    o2::tof::TOFChannelData* c = slot.getContainer();
    LOG(INFO) << "Finalize slot " << slot.getTFStart() << " <= TF <= " << slot.getTFEnd();

    // for the CCDB entry
    std::map<std::string, std::string> md;
    TimeSlewing& ts = mCalibTOFapi->getSlewParamObj(); // we take the current CCDB object, since we want to simply update the offset

#if defined(WITH_OPENMP) && !defined(__CLING__)
    if (mNThreads < 1) mNThreads = omp_get_max_threads();
    #pragma omp parallel for schedule(dynamic) num_threads(mNThreads)
#endif
    for (int sector = 0; sector < Geo::NSECTORS; sector++) {
      for (int chinsector = 0; chinsector < Geo::NPADSXSECTOR; chinsector++) {
        // make the slice of the 2D histogram so that we have the 1D of the current channel
        int ich = chinsector + sector*Geo::NPADSXSECTOR;
        auto entriesInChannel = c->integral(ich);
        if (entriesInChannel < mMinEntries) {
          LOG(DEBUG) << "channel " << ich << " will not be calibrated since it has only " << entriesInChannel << " entries (min = " << mMinEntries << ")";
          continue;
        }
        std::vector<float> fitValues;
        std::vector<float> histoValues;
        std::vector<int> entriesPerChannel = c->getEntriesPerChannel();
        if (entriesPerChannel.at(ich) == 0) {
          continue; // skip always since a channel with 0 entries is normal, it will be flagged as problematic
          if (mTest) {
            LOG(DEBUG) << "Skipping channel " << ich << " because it has zero entries, but it should not be"; // should become error!
            continue;
          } else {
            throw std::runtime_error("We found one channel with no entries, we cannot calibrate!");
          }
        }

        // more efficient way
        auto histo = c->getHisto(sector);
        for (unsigned j = chinsector; j <= chinsector; ++j) {
          for (unsigned i = 0; i < c->getNbins(); ++i) {
            const auto& v = histo.at(i, j);
            LOG(DEBUG) << "channel = " << ich << ", in sector = " << sector << " (where it is channel = " << chinsector << ") bin = " << i << " value = " << v;
            histoValues.push_back(v);
          }
        }

        double fitres = fitGaus(c->getNbins(), histoValues.data(), -(c->getRange()), c->getRange(), fitValues);

        if (fitres >= 0) {
          LOG(DEBUG) << "Channel " << ich << " :: Fit result " << fitres << " Mean = " << fitValues[1] << " Sigma = " << fitValues[2];
        } else {
          //        LOG(INFO) << "Channel " << ich << " :: Fit failed with result = " << fitres;
          continue;
        }

        if (fitValues[2] < 0) {
          fitValues[2] = -fitValues[2];
        }

        float fractionUnderPeak;
        float intmin = fitValues[1] - 5 * fitValues[2]; // mean - 5*sigma
        float intmax = fitValues[1] + 5 * fitValues[2]; // mean + 5*sigma

        if (intmin < -mRange) {
          intmin = -mRange;
        }
        if (intmax < -mRange) {
          intmax = -mRange;
        }
        if (intmin > mRange) {
          intmin = mRange;
        }
        if (intmax > mRange) {
          intmax = mRange;
        }

        fractionUnderPeak = entriesInChannel > 0 ? c->integral(ich, intmin, intmax) / entriesInChannel : 0;
        // now we need to store the results in the TimeSlewingObject
        ts.setFractionUnderPeak(ich / Geo::NPADSXSECTOR, ich % Geo::NPADSXSECTOR, fractionUnderPeak);
        ts.setSigmaPeak(ich / Geo::NPADSXSECTOR, ich % Geo::NPADSXSECTOR, abs(fitValues[2]));
        ts.updateOffsetInfo(ich, fitValues[1]);
      } // end loop channels in sector
    } // end loop over sectors
    auto clName = o2::utils::MemFileHelper::getClassName(ts);
    auto flName = o2::ccdb::CcdbApi::generateFileName(clName);
    mInfoVector.emplace_back("TOF/ChannelCalib", clName, flName, md, slot.getTFStart(), 99999999999999);
    mTimeSlewingVector.emplace_back(ts);
  }

//_____________________________________________

template class TOFChannelCalibrator<o2::dataformats::CalibInfoTOF>;
template class TOFChannelCalibrator<o2::tof::CalibInfoCluster>;


} // end namespace tof
} // end namespace o2
