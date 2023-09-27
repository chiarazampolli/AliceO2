// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GlobalTracking/MatchITSTPCQC.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "Framework/InputSpec.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "DetectorsBase/Propagator.h"
#include "SimulationDataFormat/MCUtils.h"
#include <algorithm>
#include "TGraphAsymmErrors.h"
#include "GlobalTracking/TrackCuts.h"
#include <DetectorsBase/GRPGeomHelper.h>

using namespace o2::globaltracking;
using namespace o2::mcutils;
using MCTrack = o2::MCTrackT<float>;

MatchITSTPCQC::~MatchITSTPCQC()
{

  deleteHistograms();
}

//_______________________________________________________

void MatchITSTPCQC::deleteHistograms()
{

  for (int i = 0; i < matchType::SIZE; ++i) {
    // Pt
    delete mPtNum[i];
    delete mPtDen[i];
    delete mFractionITSTPCmatch[i];
    delete mPtNum_noEta0[i];
    delete mPtDen_noEta0[i];
    delete mFractionITSTPCmatch_noEta0[i];
    delete mPtPhysPrimNum[i];
    delete mPtPhysPrimDen[i];
    delete mFractionITSTPCmatchPhysPrim[i];
    // Phi
    delete mPhiNum[i];
    delete mPhiDen[i];
    delete mFractionITSTPCmatchPhi[i];
    delete mPhiPhysPrimNum[i];
    delete mPhiPhysPrimDen[i];
    delete mFractionITSTPCmatchPhiPhysPrim[i];
    delete mPhiVsPtNum[i];
    delete mPhiVsPtDen[i];
    delete mFractionITSTPCmatchPhiVsPt[i];
    // Eta
    delete mEtaNum[i];
    delete mEtaDen[i];
    delete mFractionITSTPCmatchEta[i];
    delete mEtaPhysPrimNum[i];
    delete mEtaPhysPrimDen[i];
    delete mFractionITSTPCmatchEtaPhysPrim[i];
    delete mEtaVsPtNum[i];
    delete mEtaVsPtDen[i];
    delete mFractionITSTPCmatchEtaVsPt[i];
    // 1/Pt
    delete m1OverPtNum[i];
    delete m1OverPtDen[i];
    delete mFractionITSTPCmatch1OverPt[i];
    delete m1OverPtPhysPrimNum[i];
    delete m1OverPtPhysPrimDen[i];
    delete mFractionITSTPCmatchPhysPrim1OverPt[i];
  }

  // Residuals
  delete mResidualPt;
  delete mResidualPhi;
  delete mResidualEta;
  // Others
  delete mChi2Matching;
  delete mChi2Refit;
  delete mTimeResVsPt;
  delete mDCAr;
}

//__________________________________________________________

void MatchITSTPCQC::reset()
{
  for (int i = 0; i < matchType::SIZE; ++i) {
    // Pt
    mPtNum[i]->Reset();
    mPtDen[i]->Reset();
    mPtNum_noEta0[i]->Reset();
    mPtDen_noEta0[i]->Reset();
    mPtPhysPrimNum[i]->Reset();
    mPtPhysPrimDen[i]->Reset();
    // Phi
    mPhiNum[i]->Reset();
    mPhiDen[i]->Reset();
    mPhiPhysPrimNum[i]->Reset();
    mPhiPhysPrimDen[i]->Reset();
    mPhiVsPtNum[i]->Reset();
    mPhiVsPtDen[i]->Reset();
    // Eta
    mEtaNum[i]->Reset();
    mEtaDen[i]->Reset();
    mEtaPhysPrimNum[i]->Reset();
    mEtaPhysPrimDen[i]->Reset();
    mEtaVsPtNum[i]->Reset();
    mEtaVsPtDen[i]->Reset();
    // 1/Pt
    m1OverPtNum[i]->Reset();
    m1OverPtDen[i]->Reset();
    m1OverPtPhysPrimNum[i]->Reset();
    m1OverPtPhysPrimDen[i]->Reset();
  }
  // Residuals
  mResidualPt->Reset();
  mResidualPhi->Reset();
  mResidualEta->Reset();
  // Others
  mChi2Matching->Reset();
  mChi2Refit->Reset();
  mTimeResVsPt->Reset();
  mDCAr->Reset();
}

//__________________________________________________________
bool MatchITSTPCQC::init()
{

  std::array<std::string, 2> title{"TPC", "ITS"};
  for (int i = 0; i < matchType::SIZE; ++i) {
    // Pt
    mPtNum[i] = new TH1D(Form("mPtNum_%s", title[i].c_str()), Form("Pt distribution of ITSTPC matched tracks, wrt %s tracks; Pt [GeV/c]; dNdPt", title[i].c_str()), 100, 0.f, 20.f);
    mPtDen[i] = new TH1D(Form("mPtDen_%s", title[i].c_str()), Form("Pt distribution of %s tracks; Pt [GeV/c]; dNdPt", title[i].c_str()), 100, 0.f, 20.f);
    mFractionITSTPCmatch[i] = new TEfficiency(Form("mFractionITSTPCmatch_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks wrt %s tracks vs Pt; Pt [GeV/c]; Eff", title[i].c_str()), 100, 0.f, 20.f);
    mPtNum_noEta0[i] = new TH1D(Form("mPtNum_noEta0_%s", title[i].c_str()), Form("Pt distribution of ITSTPC matched tracks without |eta| < 0.05, wrt %s tracks; Pt [GeV/c]; dNdPt", title[i].c_str()), 100, 0.f, 20.f);
    mPtDen_noEta0[i] = new TH1D(Form("mPtDen_noEta0_%s", title[i].c_str()), Form("Pt distribution of %s tracks without |eta| < 0.05; Pt [GeV/c]; dNdPt", title[i].c_str()), 100, 0.f, 20.f);
    mFractionITSTPCmatch_noEta0[i] = new TEfficiency(Form("mFractionITSTPCmatch_noEta0_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks wrt %s tracks vs Pt without |eta| < 0.05; Pt [GeV/c]; Eff", title[i].c_str()), 100, 0.f, 20.f);
    // Phi
    mPhiNum[i] = new TH1F(Form("mPhiNum_%s", title[i].c_str()), Form("Phi distribution of ITSTPC matched tracks, wrt %s tracks; Phi [rad]; dNdPhi", title[i].c_str()), 100, 0.f, 2 * TMath::Pi());
    mPhiDen[i] = new TH1F(Form("mPhiDen_%s", title[i].c_str()), Form("Phi distribution of %s tracks; Phi [rad]; dNdPhi", title[i].c_str()), 100, 0.f, 2 * TMath::Pi());
    mFractionITSTPCmatchPhi[i] = new TEfficiency(Form("mFractionITSTPCmatchPhi_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks vs Phi; Phi [rad]; Eff", title[i].c_str()), 100, 0.f, 2 * TMath::Pi());
    mPhiVsPtNum[i] = new TH2F(Form("mPhiVsPtNum_%s", title[i].c_str()), Form("Phi vs Pt distribution of ITSTPC matched tracks; #it{p}_{T} [GeV#it{c}]; Phi [rad]; dNdPhi", title[i].c_str()), 100, 0.f, 20.f, 100, 0.f, 2 * TMath::Pi());
    mPhiVsPtDen[i] = new TH2F(Form("mPhiVsPtDen_%s", title[i].c_str()), Form("Phi vs Pt distribution of %s tracks; #it{p}_{T} [GeV#it{c}]; Phi [rad]; dNdPhi", title[i].c_str()), 100, 0.f, 20.f, 100, 0.f, 2 * TMath::Pi());
    mFractionITSTPCmatchPhiVsPt[i] = new TEfficiency(Form("mFractionITSTPCmatchPhiVsPt_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks wrt %s tracks, Phi vs Pt; #it{p}_{T} [GeV#it{c}]; Phi [rad]; Eff", title[i].c_str()), 100, 0.f, 20.f, 100, 0.f, 2 * TMath::Pi());
    // Eta
    mEtaNum[i] = new TH1F(Form("mEtaNum_%s", title[i].c_str()), Form("Eta distribution of ITSTPC matched tracks, wrt %s tracks; Eta; dNdEta", title[i].c_str()), 100, -2.f, 2.f);
    mEtaDen[i] = new TH1F(Form("mEtaDen_%s", title[i].c_str()), Form("Eta distribution of %s tracks; Eta; dNdEta", title[i].c_str()), 100, -2.f, 2.f);
    mFractionITSTPCmatchEta[i] = new TEfficiency(Form("mFractionITSTPCmatchEta_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks , wrt %s tracks, vs Eta; Eta; Eff", title[i].c_str()), 100, -2.f, 2.f);
    mEtaVsPtNum[i] = new TH2F(Form("mEtaVsPtNum_%s", title[i].c_str()), Form("Eta vs Pt distribution of ITSTPC matched tracks, wrt %s tracks; #it{p}_{T} [GeV#it{c}]; #it{p}_{T} [GeV#it{c}]; Eta", title[i].c_str()), 100, 0.f, 20.f, 100, -2.f, 2.f);
    mEtaVsPtDen[i] = new TH2F(Form("mEtaVsPtDen_%s", title[i].c_str()), Form("Eta vs Pt distribution of %s tracks; #it{p}_{T} [GeV#it{c}]; #it{p}_{T} [GeV#it{c}]; Eta", title[i].c_str()), 100, 0.f, 20.f, 100, -2.f, 2.f);
    mFractionITSTPCmatchEtaVsPt[i] = new TEfficiency(Form("mFractionITSTPCmatchEtaVsPt_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks, wrt %s tracks, Eta vs Pt; #it{p}_{T} [GeV#it{c}]; Eta; Eff", title[i].c_str()), 100, 0.f, 20.f, 100, -2.f, 2.f);
    // 1/pt
    m1OverPtNum[i] = new TH1D(Form("m1OverPtNum_%s", title[i].c_str()), Form("1/Pt distribution of matched tracks, wrt %s tracks; 1/Pt [c/GeV]; dNdPt", title[i].c_str()), 100, 0.f, 20.f);
    m1OverPtDen[i] = new TH1D(Form("m1OverPtDen_%s", title[i].c_str()), Form("1/Pt distribution of %s tracks; 1/Pt [c/GeV]; dNdPt", title[i].c_str()), 100, 0.f, 20.f);
    mFractionITSTPCmatch1OverPt[i] = new TEfficiency(Form("mFractionITSTPCmatch1OverPt_%s", title[i].c_str()), Form("Fraction of ITSTPC matched tracks vs 1/Pt, wrt %s tracks; 1/Pt [c/GeV]; Eff", title[i].c_str()), 100, 0.f, 20.f);

    // These will be empty in case of no MC info...
    mPhiPhysPrimNum[i] = new TH1F(Form("mPhiPhysPrimNum", title[i].c_str()), Form("Phi distribution of matched tracks (physical primary), wrt %s tracks; Phi [rad]; dNdPhi", title[i].c_str()), 100, 0.f, 2 * TMath::Pi());
    mPhiPhysPrimDen[i] = new TH1F(Form("mPhiPhysPrimDen", title[i].c_str()), Form("Phi distribution of %s tracks (physical primary); Phi [rad]; dNdPhi", title[i].c_str()), 100, 0.f, 2 * TMath::Pi());
    mFractionITSTPCmatchPhiPhysPrim[i] = new TEfficiency(Form("mFractionITSTPCmatchPhiPhysPrim", title[i].c_str()), Form("Fraction of ITSTPC matched tracks vs Phi (physical primary), wrt %s tracks; Phi [rad]; Eff", title[i].c_str()), 100, 0.f, 2 * TMath::Pi());
    mEtaPhysPrimNum[i] = new TH1F(Form("mEtaPhysPrimNum", title[i].c_str()), Form("Eta distribution of matched tracks (physical primary), wrt %s tracks; Eta; dNdEta", title[i].c_str()), 100, -2.f, 2.f);
    mEtaPhysPrimDen[i] = new TH1F(Form("mEtaPhysPrimDen", title[i].c_str()), Form("Eta distribution of %s tracks (physical primary); Eta; dNdEta", title[i].c_str()), 100, -2.f, 2.f);
    mFractionITSTPCmatchEtaPhysPrim[i] = new TEfficiency(Form("mFractionITSTPCmatchEtaPhysPrim", title[i].c_str()), Form("Fraction of ITSTPC matched tracks vs Eta (physical primary), wrt %s tracks; Eta; Eff", title[i].c_str()), 100, -2.f, 2.f);
  }

  mResidualPt = new TH2F("mResidualPt", "Residuals of ITS-TPC matching in #it{p}_{T}; #it{p}_{T}^{ITS-TPC} [GeV/c]; #it{p}_{T}^{ITS-TPC} - #it{p}_{T}^{TPC} [GeV/c]", 100, 0.f, 20.f, 100, -1.f, 1.f);
  mResidualPhi = new TH2F("mResidualPhi", "Residuals of ITS-TPC matching in #it{#phi}; #it{#phi}^{ITS-TPC} [rad]; #it{#phi}^{ITS-TPC} - #it{#phi}^{TPC} [rad]", 100, 0.f, 2 * TMath::Pi(), 100, -1.f, 1.f);
  mResidualEta = new TH2F("mResidualEta", "Residuals of ITS-TPC matching in #it{#eta}; #it{#eta}^{ITS-TPC}; #it{#eta}^{ITS-TPC} - #it{#eta}^{TPC}", 100, -2.f, 2.f, 100, -1.f, 1.f);
  mChi2Matching = new TH1F("mChi2Matching", "Chi2 of matching; chi2", 200, 0, 300);
  mChi2Refit = new TH1F("mChi2Refit", "Chi2 of refit; chi2", 200, 0, 300);
  mDCAr = new TH1F("mDCAr", "DCA of TPC tracks; DCAr", 200, -100, 100);

  // log binning for pT
  const Int_t nbinsPt = 100;
  const Double_t xminPt = 0.01;
  const Double_t xmaxPt = 20;
  Double_t* xbinsPt = new Double_t[nbinsPt + 1];
  Double_t xlogminPt = TMath::Log10(xminPt);
  Double_t xlogmaxPt = TMath::Log10(xmaxPt);
  Double_t dlogxPt = (xlogmaxPt - xlogminPt) / nbinsPt;
  for (int i = 0; i <= nbinsPt; i++) {
    Double_t xlogPt = xlogminPt + i * dlogxPt;
    xbinsPt[i] = TMath::Exp(TMath::Log(10) * xlogPt);
  }
  mTimeResVsPt = new TH2F("mTimeResVsPt", "Time resolution vs Pt; Pt [GeV/c]; time res [us]", nbinsPt, xbinsPt, 100, 0.f, 2.f);

  for (int i = 0; i < matchType::SIZE; ++i) {
    mPtPhysPrimNum[i] = new TH1F(Form("mPtPhysPrimNum", title[i].c_str()), Form("Pt distribution of matched tracks (physical primary), wrt %s tracks; Pt [GeV/c]; dNdPt", title[i].c_str()), nbinsPt, xbinsPt);
    mPtPhysPrimDen[i] = new TH1F(Form("mPtPhysPrimDen", title[i].c_str()), Form("Pt distribution of %s tracks (physical primary); Pt [GeV/c]; dNdPt", title[i].c_str()), nbinsPt, xbinsPt);
    mFractionITSTPCmatchPhysPrim[i] = new TEfficiency(Form("mFractionITSTPCmatchPhysPrim", title[i].c_str()), Form("Fraction of ITSTPC matched tracks vs Pt (physical primary), wrt %s tracks; Pt [GeV/c]; Eff", title[i].c_str()), nbinsPt, xbinsPt);
    m1OverPtPhysPrimNum[i] = new TH1D(Form("m1OverPtPhysPrimNum", title[i].c_str()), Form("1/Pt distribution of matched tracks (physical primary), wrt %s tracks; 1/Pt [c/GeV]; dNd1/Pt", title[i].c_str()), nbinsPt, xbinsPt);
    m1OverPtPhysPrimDen[i] = new TH1D(Form("m1OverPtPhysPrimDen", title[i].c_str()), Form("1/PtPt distribution of %s tracks (physical primary); 1/Pt [c/GeV]; dNd1/Pt", title[i].c_str()), nbinsPt, xbinsPt);
    mFractionITSTPCmatchPhysPrim1OverPt[i] = new TEfficiency(Form("mFractionITSTPCmatchPhysPrim1OverPt", title[i].c_str()), Form("Fraction of ITSTPC matched tracks vs 1/Pt (physical primary), wrt %s tracks; 1/Pt [c/GeV]; Eff", title[i].c_str()), nbinsPt, xbinsPt);

    // some extra settings
    mPtNum[i]->Sumw2();
    mPtDen[i]->Sumw2();
    mPtNum_noEta0[i]->Sumw2();
    mPtDen_noEta0[i]->Sumw2();
    mPhiNum[i]->Sumw2();
    mPhiDen[i]->Sumw2();
    mPhiVsPtNum[i]->Sumw2();
    mPhiVsPtDen[i]->Sumw2();
    mPtPhysPrimNum[i]->Sumw2();
    mPtPhysPrimDen[i]->Sumw2();
    mPhiPhysPrimNum[i]->Sumw2();
    mPhiPhysPrimDen[i]->Sumw2();
    mEtaNum[i]->Sumw2();
    mEtaDen[i]->Sumw2();
    mEtaPhysPrimNum[i]->Sumw2();
    mEtaPhysPrimDen[i]->Sumw2();
    mEtaVsPtNum[i]->Sumw2();
    mEtaVsPtDen[i]->Sumw2();

    m1OverPtNum[i]->Sumw2();
    m1OverPtDen[i]->Sumw2();
    m1OverPtPhysPrimNum[i]->Sumw2();
    m1OverPtPhysPrimDen[i]->Sumw2();

    mPtNum_noEta0[i]->SetOption("logy");
    mPtDen_noEta0[i]->SetOption("logy");
    mPtNum[i]->SetOption("logy");
    mPtDen[i]->SetOption("logy");

    mPtNum[i]->GetYaxis()->SetTitleOffset(1.4);
    mPtDen[i]->GetYaxis()->SetTitleOffset(1.4);
    mPtNum_noEta0[i]->GetYaxis()->SetTitleOffset(1.4);
    mPtDen_noEta0[i]->GetYaxis()->SetTitleOffset(1.4);
    mEtaNum[i]->GetYaxis()->SetTitleOffset(1.4);
    mEtaDen[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  mChi2Matching->SetOption("logy");
  mChi2Refit->SetOption("logy");
  mTimeResVsPt->SetOption("colz logz logy logx");

  mChi2Matching->GetYaxis()->SetTitleOffset(1.4);
  mChi2Refit->GetYaxis()->SetTitleOffset(1.4);
  mTimeResVsPt->GetYaxis()->SetTitleOffset(1.4);

  if (mUseMC) {
    mcReader.initFromDigitContext("collisioncontext.root");
  }

  return true;
}

//__________________________________________________________

void MatchITSTPCQC::initDataRequest()
{

  // initialize data request, if it was not already done

  mSrc &= mAllowedSources;

  if (mSrc[GID::Source::ITSTPC] == 0 || mSrc[GID::Source::TPC] == 0 || mSrc[GID::Source::ITS] == 0) {
    LOG(fatal) << "We cannot do ITSTPC QC, some sources are missing, check sources in " << mSrc;
  }

  mDataRequest = std::make_shared<o2::globaltracking::DataRequest>();
  mDataRequest->requestTracks(mSrc, mUseMC);
}

//__________________________________________________________

void MatchITSTPCQC::run(o2::framework::ProcessingContext& ctx)
{

  // Getting the B field
  mBz = o2::base::Propagator::Instance()->getNominalBz();

  static int evCount = 0;
  mRecoCont.collectData(ctx, *mDataRequest.get());
  mTPCTracks = mRecoCont.getTPCTracks();
  mITSTracks = mRecoCont.getITSTracks();
  mITSTPCTracks = mRecoCont.getTPCITSTracks();

  LOG(debug) << "****** Number of found ITSTPC tracks = " << mITSTPCTracks.size();
  LOG(debug) << "****** Number of found TPC    tracks = " << mTPCTracks.size();
  LOG(debug) << "****** Number of found ITS    tracks = " << mITSTracks.size();

  // cache selection for TPC and ITS tracks
  std::vector<bool> isTPCTrackSelectedEntry(mTPCTracks.size(), false);
  std::vector<bool> isITSTrackSelectedEntry(mITSTracks.size(), false);
  TrackCuts cuts;

  for (size_t itrk = 0; itrk < mTPCTracks.size(); ++itrk) {
    auto const& trkTpc = mTPCTracks[itrk];
    // if (selectTrack(trkTpc)) {
    //   isTPCTrackSelectedEntry[itrk] = true;
    // }
    o2::dataformats::GlobalTrackID id(itrk, GID::TPC);
    if (cuts.isSelected(id, mRecoCont)) {
      isTPCTrackSelectedEntry[itrk] = true;
    }
  }

  for (size_t itrk = 0; itrk < mITSTracks.size(); ++itrk) {
    auto const& trkIts = mITSTracks[itrk];
    o2::dataformats::GlobalTrackID id(itrk, GID::ITS);
    if (cuts.isSelected(id, mRecoCont)) {
      isITSTrackSelectedEntry[itrk] = true;
    }
  }

  // numerator + eta, chi2...
  if (mUseMC) {
    for (int i = 0; i < matchType::SIZE; ++i) {
      mMapLabels[i].clear();
    }
    for (int itrk = 0; itrk < static_cast<int>(mITSTPCTracks.size()); ++itrk) {
      auto const& trk = mITSTPCTracks[itrk];
      auto idxTrkTpc = trk.getRefTPC().getIndex();
      if (isTPCTrackSelectedEntry[idxTrkTpc] == true) {
        auto lbl = mRecoCont.getTrackMCLabel({(unsigned int)(itrk), GID::Source::ITSTPC});
        if (mMapLabels[matchType::TPC].find(lbl) == mMapLabels[matchType::TPC].end()) {
          int source = lbl.getSourceID();
          int event = lbl.getEventID();
          const std::vector<o2::MCTrack>& pcontainer = mcReader.getTracks(source, event);
          const o2::MCTrack& p = pcontainer[lbl.getTrackID()];
          if (MCTrackNavigator::isPhysicalPrimary(p, pcontainer)) {
            mMapLabels[matchType::TPC].insert({lbl, {itrk, true}});
          } else {
            mMapLabels[matchType::TPC].insert({lbl, {itrk, false}});
          }
        } else {
          // winner (if more tracks have the same label) has the highest pt
          if (mITSTPCTracks[mMapLabels[matchType::TPC].at(lbl).mIdx].getPt() < trk.getPt()) {
            mMapLabels[matchType::TPC].at(lbl).mIdx = itrk;
          }
        }
      }
      auto idxTrkIts = trk.getRefITS().getIndex();
      /*
  if (trk.getRefITS().getSource() == GID::Source::ITSAB) {
  // track from AfterBurner, skipping
  continue;
      }
      */
      if (isITSTrackSelectedEntry[idxTrkIts] == true) {
        auto lbl = mRecoCont.getTrackMCLabel({(unsigned int)(itrk), GID::Source::ITSTPC});
        if (mMapLabels[matchType::ITS].find(lbl) == mMapLabels[matchType::ITS].end()) {
          int source = lbl.getSourceID();
          int event = lbl.getEventID();
          const std::vector<o2::MCTrack>& pcontainer = mcReader.getTracks(source, event);
          const o2::MCTrack& p = pcontainer[lbl.getTrackID()];
          if (MCTrackNavigator::isPhysicalPrimary(p, pcontainer)) {
            mMapLabels[matchType::ITS].insert({lbl, {itrk, true}});
          } else {
            mMapLabels[matchType::ITS].insert({lbl, {itrk, false}});
          }
        } else {
          // winner (if more tracks have the same label) has the highest pt
          if (mITSTPCTracks[mMapLabels[matchType::ITS].at(lbl).mIdx].getPt() < trk.getPt()) {
            mMapLabels[matchType::ITS].at(lbl).mIdx = itrk;
          }
        }
      }
    }
    LOG(info) << "number of entries in map for nominator (without duplicates) = " << mMapLabels.size();
    // now we use only the tracks in the map to fill the histograms (--> tracks have passed the
    // track selection and there are no duplicated tracks wrt the same MC label)
    for (int i = 0; i < matchType::SIZE; ++i) {
      for (auto const& el : mMapLabels[i]) {
        auto const& trk = mITSTPCTracks[el.second.mIdx];
        o2::track::TrackParCov trkDen;
        if (i == matchType::TPC) {
          trkDen = mTPCTracks[trk.getRefTPC()];
        } else {
          trkDen = mITSTracks[trk.getRefITS()];
        }
        mPtNum[i]->Fill(trkDen.getPt());
        if (std::abs(trkDen.getEta()) > 0.05) {
          mPtNum_noEta0[i]->Fill(trkDen.getPt());
        }
        mPhiNum[i]->Fill(trkDen.getPhi());
        mPhiVsPtNum[i]->Fill(trkDen.getPt(), trkDen.getPhi());
        mEtaNum[i]->Fill(trkDen.getEta());
        mEtaVsPtNum[i]->Fill(trkDen.getPt(), trkDen.getEta());
        m1OverPtNum[i]->Fill(trkDen.getPtInv());
        // we fill also the denominator
        mPtDen[i]->Fill(trkDen.getPt());
        if (std::abs(trkDen.getEta()) > 0.05) {
          mPtDen_noEta0[i]->Fill(trkDen.getPt());
        }
        mPhiDen[i]->Fill(trkDen.getPhi());
        mPhiVsPtDen[i]->Fill(trkDen.getPt(), trkDen.getPhi());
        mEtaDen[i]->Fill(trkDen.getEta());
        mEtaVsPtDen[i]->Fill(trkDen.getPt(), trkDen.getEta());
        m1OverPtDen[i]->Fill(trkDen.getPtInv());
        if (el.second.mIsPhysicalPrimary) {
          mPtPhysPrimNum[i]->Fill(trkDen.getPt());
          mPhiPhysPrimNum[i]->Fill(trkDen.getPhi());
          mEtaPhysPrimNum[i]->Fill(trkDen.getEta());
          m1OverPtPhysPrimNum[i]->Fill(trkDen.getPtInv());
          // we fill also the denominator
          mPtPhysPrimDen[i]->Fill(trkDen.getPt());
          mPhiPhysPrimDen[i]->Fill(trkDen.getPhi());
          mEtaPhysPrimDen[i]->Fill(trkDen.getEta());
          m1OverPtPhysPrimDen[i]->Fill(trkDen.getPtInv());
        }
        ++mNITSTPCSelectedTracks[i];
      }
    }
  }
  int iITSTPC = 0;
  for (auto const& trk : mITSTPCTracks) {
    if (trk.getRefTPC().getIndex() >= mTPCTracks.size()) {
      LOG(fatal) << "******************** ATTENTION! for TPC track associated to matched track: idx = " << trk.getRefTPC().getIndex() << ", size of container = " << mTPCTracks.size() << " in TF " << evCount;
    }
    if (trk.getRefITS().getIndex() >= mITSTracks.size()) {
      LOG(fatal) << "******************** ATTENTION! for ITS track associated to matched track: idx = " << trk.getRefITS().getIndex() << ", size of container = " << mITSTracks.size() << " in TF " << evCount;
    }
    std::array<std::string, 2> title{"TPC", "ITS"};
    for (int i = 0; i < matchType::SIZE; ++i) {
      o2::track::TrackParCov trkRef;
      int idxTrkRef;
      bool fillHisto = false;
      if (i == matchType::TPC) {
        trkRef = mTPCTracks[trk.getRefTPC()];
        idxTrkRef = trk.getRefTPC().getIndex();
        if (isTPCTrackSelectedEntry[idxTrkRef] == true) {
          fillHisto = true;
          ++mNITSTPCSelectedTracks[i];
        }
      } else {
        trkRef = mITSTracks[trk.getRefITS()];
        idxTrkRef = trk.getRefITS().getIndex();
        LOG(debug) << "Checking track (ITS) with id " << idxTrkRef << " for ITSTPC track " << iITSTPC << " and pt = " << trkRef.getPt();
        if (trk.getRefITS().getSource() == GID::ITSAB) {
          // do not use afterburner tracks
          LOG(debug) << "Track (ITS) with id " << idxTrkRef << " for ITSTPC track " << iITSTPC << " is from afterburner";
          continue;
        }
        if (isITSTrackSelectedEntry[idxTrkRef] == true) {
          LOG(debug) << "Track was selected (ITS), with id " << idxTrkRef << " for ITSTPC track " << iITSTPC << " , we keep it in the numerator, pt = " << trkRef.getPt();
          fillHisto = true;
          ++mNITSTPCSelectedTracks[i];
        } else {
          LOG(debug) << "Track was not selected (ITS), with id " << idxTrkRef << " for ITSTPC track " << iITSTPC << " , we don't keep it in the numerator, pt = " << trkRef.getPt();
        }
      }
      if (fillHisto == true) {
        if (!mUseMC) {
          LOG(debug) << "Filling num (" << title[i] << ") with track with id " << idxTrkRef << " for ITSTPC track " << iITSTPC << " with pt = " << trkRef.getPt();
          mPtNum[i]->Fill(trkRef.getPt());
          if (std::abs(trkRef.getEta()) > 0.05) {
            mPtNum_noEta0[i]->Fill(trkRef.getPt());
          }
          mPhiNum[i]->Fill(trkRef.getPhi());
          mPhiVsPtNum[i]->Fill(trkRef.getPt(), trkRef.getPhi());
          mEtaNum[i]->Fill(trkRef.getEta());
          mEtaVsPtNum[i]->Fill(trkRef.getPt(), trkRef.getEta());
          m1OverPtNum[i]->Fill(trkRef.getPtInv());
        }
        if (i == matchType::TPC) {
          mResidualPt->Fill(trk.getPt(), trk.getPt() - trkRef.getPt());
          mResidualPhi->Fill(trk.getPhi(), trk.getPhi() - trkRef.getPhi());
          mResidualEta->Fill(trk.getEta(), trk.getEta() - trkRef.getEta());
          mChi2Matching->Fill(trk.getChi2Match());
          mChi2Refit->Fill(trk.getChi2Refit());
          mTimeResVsPt->Fill(trkRef.getPt(), trk.getTimeMUS().getTimeStampError());
          math_utils::Point3D<float> v{};
          std::array<float, 2> dca;
          if (trkRef.propagateParamToDCA(v, mBz, &dca)) {
            mDCAr->Fill(dca[0]);
          }
          LOG(debug) << "*** chi2Matching = " << trk.getChi2Match() << ", chi2refit = " << trk.getChi2Refit() << ", timeResolution = " << trk.getTimeMUS().getTimeStampError();
        }
      } else {
        LOG(debug) << "Not filling num (" << title[i] << ") for ITSTPC track " << iITSTPC << " for track with pt " << trkRef.getPt();
      }
    }
    ++iITSTPC;
  }

  // now filling the denominator for the efficiency calculation
  if (mUseMC) {
    for (int i = 0; i < matchType::SIZE; ++i) {
      mMapRefLabels[i].clear();
    }
    // filling the map where we store for each MC label, the track id of the reconstructed
    // track with the highest number of TPC clusters
    for (int itrk = 0; itrk < static_cast<int>(mTPCTracks.size()); ++itrk) {
      auto const& trk = mTPCTracks[itrk];
      if (isTPCTrackSelectedEntry[itrk] == true) {
        auto lbl = mRecoCont.getTrackMCLabel({(unsigned int)(itrk), GID::Source::TPC});
        if (mMapLabels[matchType::TPC].find(lbl) != mMapLabels[matchType::TPC].end()) {
          // the track was already added to the denominator
          continue;
        }
        if (mMapRefLabels[matchType::TPC].find(lbl) == mMapRefLabels[matchType::TPC].end()) {
          int source = lbl.getSourceID();
          int event = lbl.getEventID();
          const std::vector<o2::MCTrack>& pcontainer = mcReader.getTracks(source, event);
          const o2::MCTrack& p = pcontainer[lbl.getTrackID()];
          if (MCTrackNavigator::isPhysicalPrimary(p, pcontainer)) {
            mMapRefLabels[matchType::TPC].insert({lbl, {itrk, true}});
          } else {
            mMapRefLabels[matchType::TPC].insert({lbl, {itrk, false}});
          }
        } else {
          // winner (if more tracks have the same label) has the highest number of TPC clusters
          if (mTPCTracks[mMapRefLabels[matchType::TPC].at(lbl).mIdx].getNClusters() < trk.getNClusters()) {
            mMapRefLabels[matchType::TPC].at(lbl).mIdx = itrk;
          }
        }
      }
    }
    // same for ITS
    // filling the map where we store for each MC label, the track id of the reconstructed
    // track with the highest number of ITS clusters
    for (int itrk = 0; itrk < static_cast<int>(mITSTracks.size()); ++itrk) {
      auto const& trk = mITSTracks[itrk];
      if (isITSTrackSelectedEntry[itrk] == true) {
        auto lbl = mRecoCont.getTrackMCLabel({(unsigned int)(itrk), GID::Source::ITS});
        if (mMapLabels[matchType::ITS].find(lbl) != mMapLabels[matchType::ITS].end()) {
          // the track was already added to the denominator
          continue;
        }
        if (mMapRefLabels[matchType::ITS].find(lbl) == mMapRefLabels[matchType::ITS].end()) {
          int source = lbl.getSourceID();
          int event = lbl.getEventID();
          const std::vector<o2::MCTrack>& pcontainer = mcReader.getTracks(source, event);
          const o2::MCTrack& p = pcontainer[lbl.getTrackID()];
          if (MCTrackNavigator::isPhysicalPrimary(p, pcontainer)) {
            mMapRefLabels[matchType::ITS].insert({lbl, {itrk, true}});
          } else {
            mMapRefLabels[matchType::ITS].insert({lbl, {itrk, false}});
          }
        } else {
          // winner (if more tracks have the same label) has the highest number of ITS clusters
          if (mITSTracks[mMapRefLabels[matchType::ITS].at(lbl).mIdx].getNClusters() < trk.getNClusters()) {
            mMapRefLabels[matchType::ITS].at(lbl).mIdx = itrk;
          }
        }
      }
    }
    LOG(info) << "number of entries in map for denominator of TPC tracks (without duplicates) = " << mMapRefLabels[matchType::TPC].size() + mMapLabels[matchType::TPC].size();
    LOG(info) << "number of entries in map for denominator of ITS tracks (without duplicates) = " << mMapRefLabels[matchType::ITS].size() + mMapLabels[matchType::ITS].size();
    // now we use only the tracks in the map to fill the histograms (--> tracks have passed the
    // track selection and there are no duplicated tracks wrt the same MC label)
    for (auto const& el : mMapRefLabels[matchType::TPC]) {
      auto const& trk = mTPCTracks[el.second.mIdx];
      mPtDen[matchType::TPC]->Fill(trk.getPt());
      if (std::abs(trk.getEta()) > 0.05) {
        mPtDen_noEta0[matchType::TPC]->Fill(trk.getPt());
      }
      mPhiDen[matchType::TPC]->Fill(trk.getPhi());
      mPhiVsPtDen[matchType::TPC]->Fill(trk.getPt(), trk.getPhi());
      mEtaDen[matchType::TPC]->Fill(trk.getEta());
      mEtaVsPtDen[matchType::TPC]->Fill(trk.getPt(), trk.getEta());
      m1OverPtDen[matchType::TPC]->Fill(trk.getPtInv());
      if (el.second.mIsPhysicalPrimary) {
        mPtPhysPrimDen[matchType::TPC]->Fill(trk.getPt());
        mPhiPhysPrimDen[matchType::TPC]->Fill(trk.getPhi());
        mEtaPhysPrimDen[matchType::TPC]->Fill(trk.getEta());
        m1OverPtPhysPrimDen[matchType::TPC]->Fill(trk.getPtInv());
      }
      ++mNTPCSelectedTracks;
    }
    for (auto const& el : mMapRefLabels[matchType::ITS]) {
      auto const& trk = mITSTracks[el.second.mIdx];
      mPtDen[matchType::ITS]->Fill(trk.getPt());
      if (std::abs(trk.getEta()) > 0.05) {
        mPtDen_noEta0[matchType::ITS]->Fill(trk.getPt());
      }
      mPhiDen[matchType::ITS]->Fill(trk.getPhi());
      mPhiVsPtDen[matchType::ITS]->Fill(trk.getPt(), trk.getPhi());
      mEtaDen[matchType::ITS]->Fill(trk.getEta());
      mEtaVsPtDen[matchType::ITS]->Fill(trk.getPt(), trk.getEta());
      m1OverPtDen[matchType::ITS]->Fill(trk.getPtInv());
      if (el.second.mIsPhysicalPrimary) {
        mPtPhysPrimDen[matchType::ITS]->Fill(trk.getPt());
        mPhiPhysPrimDen[matchType::ITS]->Fill(trk.getPhi());
        mEtaPhysPrimDen[matchType::ITS]->Fill(trk.getEta());
        m1OverPtPhysPrimDen[matchType::ITS]->Fill(trk.getPtInv());
      }
      ++mNITSSelectedTracks;
    }
  } else {
    // if we are in data, we loop over all tracks (no check on the label)
    for (size_t itrk = 0; itrk < mTPCTracks.size(); ++itrk) {
      auto const& trk = mTPCTracks[itrk];
      if (isTPCTrackSelectedEntry[itrk] == true) {
        LOG(debug) << "Filling den (TPC) with track with pt = " << trk.getPt();
        mPtDen[matchType::TPC]->Fill(trk.getPt());
        if (std::abs(trk.getEta()) > 0.05) {
          mPtDen_noEta0[matchType::TPC]->Fill(trk.getPt());
        }
        mPhiDen[matchType::TPC]->Fill(trk.getPhi());
        mPhiVsPtDen[matchType::TPC]->Fill(trk.getPt(), trk.getPhi());
        mEtaDen[matchType::TPC]->Fill(trk.getEta());
        mEtaVsPtDen[matchType::TPC]->Fill(trk.getPt(), trk.getEta());
        m1OverPtDen[matchType::TPC]->Fill(trk.getPtInv());
        ++mNTPCSelectedTracks;
      }
    }
    for (size_t itrk = 0; itrk < mITSTracks.size(); ++itrk) {
      auto const& trk = mITSTracks[itrk];
      LOG(debug) << "Checking den for track (ITS) " << itrk << " with pt " << trk.getPt();
      if (isITSTrackSelectedEntry[itrk] == true) {
        LOG(debug) << "Filling den for track (ITS) " << itrk << " with pt = " << trk.getPt();
        mPtDen[matchType::ITS]->Fill(trk.getPt());
        if (std::abs(trk.getEta()) > 0.05) {
          mPtDen_noEta0[matchType::ITS]->Fill(trk.getPt());
        }
        mPhiDen[matchType::ITS]->Fill(trk.getPhi());
        mPhiVsPtDen[matchType::ITS]->Fill(trk.getPt(), trk.getPhi());
        mEtaDen[matchType::ITS]->Fill(trk.getEta());
        mEtaVsPtDen[matchType::ITS]->Fill(trk.getPt(), trk.getEta());
        m1OverPtDen[matchType::ITS]->Fill(trk.getPtInv());
        ++mNITSSelectedTracks;
      } else {
        LOG(debug) << "Not filling for this track (ITS) " << itrk << " with pt = " << trk.getPt();
      }
    }
  }
  evCount++;
}

//__________________________________________________________

bool MatchITSTPCQC::selectTrack(o2::tpc::TrackTPC const& track)
{

  if (track.getPt() < mPtCut) {
    return false;
  }
  if (std::abs(track.getEta()) > mEtaCut) {
    return false;
  }
  if (track.getNClusters() < mNTPCClustersCut) {
    return false;
  }

  math_utils::Point3D<float> v{};
  std::array<float, 2> dca;
  if (!(const_cast<o2::tpc::TrackTPC&>(track).propagateParamToDCA(v, mBz, &dca, mDCACut)) || std::abs(dca[0]) > mDCACutY) {
    return false;
  }

  return true;
}

//__________________________________________________________

void MatchITSTPCQC::finalize()
{

  std::array<std::string, 2> title{"TPC", "ITS"};

  // first we use denominators and nominators to set the TEfficiency; later they are scaled

  // some checks
  for (int i = 0; i < matchType::SIZE; ++i) {
    for (int i = 0; i < mPtDen[i]->GetNbinsX(); ++i) {
      if (mPtDen[i]->GetBinContent(i + 1) < mPtNum[i]->GetBinContent(i + 1)) {
        LOG(error) << title[i] << ": bin " << i + 1 << " in [" << mPtNum[i]->GetBinLowEdge(i + 1) << " , " << mPtNum[i]->GetBinLowEdge(i + 1) + mPtNum[i]->GetBinWidth(i + 1) << "]: mPtDen[i] = " << mPtDen[i]->GetBinContent(i + 1) << ", mPtNum[i] = " << mPtNum[i]->GetBinContent(i + 1);
      }
    }
    for (int i = 0; i < mPtDen_noEta0[i]->GetNbinsX(); ++i) {
      if (mPtDen_noEta0[i]->GetBinContent(i + 1) < mPtNum_noEta0[i]->GetBinContent(i + 1)) {
        LOG(error) << title[i] << ": bin " << i + 1 << " in [" << mPtNum_noEta0[i]->GetBinLowEdge(i + 1) << " , " << mPtNum_noEta0[i]->GetBinLowEdge(i + 1) + mPtNum_noEta0[i]->GetBinWidth(i + 1) << "]: mPtDen_noEta0[i] = " << mPtDen_noEta0[i]->GetBinContent(i + 1) << ", mPtNum_noEta0[i] = " << mPtNum_noEta0[i]->GetBinContent(i + 1);
      }
    }
    for (int i = 0; i < mPhiDen[i]->GetNbinsX(); ++i) {
      if (mPhiDen[i]->GetBinContent(i + 1) < mPhiNum[i]->GetBinContent(i + 1)) {
        LOG(error) << title[i] << ": bin " << i + 1 << " in [" << mPhiNum[i]->GetBinLowEdge(i + 1) << " , " << mPhiNum[i]->GetBinLowEdge(i + 1) + mPhiNum[i]->GetBinWidth(i + 1) << "]: mPhiDen[i] = " << mPhiDen[i]->GetBinContent(i + 1) << ", mPhiNum[i] = " << mPhiNum[i]->GetBinContent(i + 1);
      }
    }
    for (int i = 0; i < mEtaDen[i]->GetNbinsX(); ++i) {
      if (mEtaDen[i]->GetBinContent(i + 1) < mEtaNum[i]->GetBinContent(i + 1)) {
        LOG(error) << title[i] << ": bin " << i + 1 << " in [" << mEtaNum[i]->GetBinLowEdge(i + 1) << " , " << mEtaNum[i]->GetBinLowEdge(i + 1) + mEtaNum[i]->GetBinWidth(i + 1) << "]: mEtaDen[i] = " << mEtaDen[i]->GetBinContent(i + 1) << ", mEtaNum[i] = " << mEtaNum[i]->GetBinContent(i + 1);
      }
    }

    // filling the efficiency
    setEfficiency(mFractionITSTPCmatch[i], mPtNum[i], mPtDen[i]);
    setEfficiency(mFractionITSTPCmatch_noEta0[i], mPtNum_noEta0[i], mPtDen_noEta0[i]);
    setEfficiency(mFractionITSTPCmatchPhi[i], mPhiNum[i], mPhiDen[i]);
    setEfficiency(mFractionITSTPCmatchEta[i], mEtaNum[i], mEtaDen[i]);
    setEfficiency(mFractionITSTPCmatchPhiVsPt[i], mPhiVsPtNum[i], mPhiVsPtDen[i], true);
    setEfficiency(mFractionITSTPCmatchEtaVsPt[i], mEtaVsPtNum[i], mEtaVsPtDen[i], true);
    setEfficiency(mFractionITSTPCmatch1OverPt[i], m1OverPtNum[i], m1OverPtDen[i]);
    if (mUseMC) {
      setEfficiency(mFractionITSTPCmatchPhysPrim[i], mPtPhysPrimNum[i], mPtPhysPrimDen[i]);
      setEfficiency(mFractionITSTPCmatchPhiPhysPrim[i], mPhiPhysPrimNum[i], mPhiPhysPrimDen[i]);
      setEfficiency(mFractionITSTPCmatchEtaPhysPrim[i], mEtaPhysPrimNum[i], mEtaPhysPrimDen[i]);
      setEfficiency(mFractionITSTPCmatchPhysPrim1OverPt[i], m1OverPtPhysPrimNum[i], m1OverPtPhysPrimDen[i]);
    }
  }
  /*
  mPtTPC->Scale(scaleFactTPC);
  mPt->Scale(scaleFactITSTPC);
  mPhiTPC->Scale(scaleFactTPC);
  mPhi->Scale(scaleFactITSTPC);
  if (mUseMC) {
    mPtTPCPhysPrim->Scale(scaleFactTPC);
    mPtPhysPrim->Scale(scaleFactITSTPC);
    mPhiTPCPhysPrim->Scale(scaleFactTPC);
    mPhiPhysPrim->Scale(scaleFactITSTPC);
  }
  mEta->Scale(scaleFactITSTPC);
  mChi2Matching->Scale(scaleFactITSTPC);
  mChi2Refit->Scale(scaleFactITSTPC);
  //mTimeResVsPt->Scale(scaleFactITSTPC); // if to few entries, one sees nothing after normalization --> let's not normalize
  */
}

//__________________________________________________________

void MatchITSTPCQC::setEfficiency(TEfficiency* eff, TH1* hnum, TH1* hden, bool is2D)
{
  if (!eff) {
    LOG(fatal) << "Cannot get TEfficiency object ";
  }
  if (!hnum) {
    LOG(fatal) << "Cannot get numerator histogram for TEfficiency object " << eff->GetName();
  }
  if (!hden) {
    LOG(fatal) << "Cannot get denominator histogram for TEfficiency object " << eff->GetName();
  }
  LOG(info) << "Setting efficiency " << eff->GetName() << " from " << hnum->GetName() << " and " << hden->GetName();

  // we need to force to replace the total histogram, otherwise it will compare it to the previous passed one, and it might get an error of inconsistency in the bin contents
  if constexpr (0) { // checking
    LOG(info) << "Num " << hnum->GetName() << " " << hnum->GetNbinsX() << " " << hnum->GetNbinsY();
    LOG(info) << "Den " << hden->GetName() << " " << hden->GetNbinsX() << " " << hden->GetNbinsY();
    for (int i = 1; i <= hden->GetNbinsX(); i++) {
      if (hden->GetBinContent(i) < hnum->GetBinContent(i)) {
        LOG(warning) << "bin " << i << " den: " << hden->GetBinContent(i) << " < num: " << hnum->GetBinContent(i) << " should be the opposite";
      }
    }
  }
  // we need to force to replace the total histogram, otherwise it will compare it to the previous passed one, and it might get an error of inconsistency in the bin contents
  if (!eff->SetTotalHistogram(*hden, "f")) {
    LOG(fatal) << "Something went wrong when defining the efficiency denominator " << eff->GetName() << " from " << hnum->GetName();
  }
  if (!eff->SetPassedHistogram(*hnum, "")) {
    LOG(fatal) << "Something went wrong when defining the efficiency numerator " << eff->GetName() << " from " << hnum->GetName();
  }
  if (is2D) {
    eff->SetTitle(Form("%s;%s;%s;%s", eff->GetTitle(), hnum->GetXaxis()->GetTitle(), hnum->GetYaxis()->GetTitle(), "Efficiency"));
  } else {
    eff->SetTitle(Form("%s;%s;%s", eff->GetTitle(), hnum->GetXaxis()->GetTitle(), "Efficiency"));
  }
}

//__________________________________________________________

void MatchITSTPCQC::getHistos(TObjArray& objar)
{

  for (int i = 0; i < matchType::SIZE; ++i) {
    objar.Add(mPtNum[i]);
    objar.Add(mPtDen[i]);
    objar.Add(mFractionITSTPCmatch[i]);

    objar.Add(mPtNum_noEta0[i]);
    objar.Add(mPtDen_noEta0[i]);
    objar.Add(mFractionITSTPCmatch_noEta0[i]);

    objar.Add(mPtPhysPrimNum[i]);
    objar.Add(mPtPhysPrimDen[i]);
    objar.Add(mFractionITSTPCmatchPhysPrim[i]);

    objar.Add(mPhiNum[i]);
    objar.Add(mPhiDen[i]);
    objar.Add(mFractionITSTPCmatchPhi[i]);

    objar.Add(mPhiPhysPrimNum[i]);
    objar.Add(mPhiPhysPrimDen[i]);
    objar.Add(mFractionITSTPCmatchPhiPhysPrim[i]);

    objar.Add(mPhiVsPtNum[i]);
    objar.Add(mPhiVsPtDen[i]);
    objar.Add(mFractionITSTPCmatchPhiVsPt[i]);

    objar.Add(mEtaNum[i]);
    objar.Add(mEtaDen[i]);
    objar.Add(mFractionITSTPCmatchEta[i]);

    objar.Add(mEtaPhysPrimNum[i]);
    objar.Add(mEtaPhysPrimDen[i]);
    objar.Add(mFractionITSTPCmatchEtaPhysPrim[i]);

    objar.Add(mEtaVsPtNum[i]);
    objar.Add(mEtaVsPtDen[i]);
    objar.Add(mFractionITSTPCmatchEtaVsPt[i]);

    objar.Add(m1OverPtNum[i]);
    objar.Add(m1OverPtDen[i]);
    objar.Add(mFractionITSTPCmatch1OverPt[i]);

    objar.Add(m1OverPtPhysPrimNum[i]);
    objar.Add(m1OverPtPhysPrimDen[i]);
    objar.Add(mFractionITSTPCmatchPhysPrim1OverPt[i]);
  }
  objar.Add(mChi2Matching);
  objar.Add(mChi2Refit);
  objar.Add(mTimeResVsPt);
  objar.Add(mResidualPt);
  objar.Add(mResidualPhi);
  objar.Add(mResidualEta);
  objar.Add(mDCAr);
}
