//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TH1D.h"
#include "TNamed.h"

//Local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/returnRootFileContentsList.h"

//Define PI
#define PI 3.14159265359

int DijetImbalanceRatio(const std::string flatPthatFileName, const std::string stagPthatFileName)
{
  bool goodFlatFile = checkFile(flatPthatFileName);
  bool goodStagFile = checkFile(stagPthatFileName);

  if (!goodFlatFile || !goodStagFile) {
    if (!goodFlatFile)
      std::cout << "Input flat pTHat file \'" << flatPthatFileName << "\' is invalid. return 1" << std::endl;
    if (!goodStagFile)
      std::cout << "Input stag pTHat file \'" << stagPthatFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime *date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");

  const std::string outFileName = "output/statComp_" + dateStr + ".root";
  TFile *outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  TFile *inStagFile_p = new TFile(stagPthatFileName.c_str(), "READ");
  std::vector < std::string > tnamedPthat = returnRootFileContentsList(inStagFile_p, "TNamed", "crossSecti");

  //Grabbing cross section values for weighting internally
  const Int_t nPthatFiles = tnamedPthat.size();
  Double_t pthatFiles[nPthatFiles + 1];
  Double_t xSections[nPthatFiles + 1];
  for (int pI = 0; pI < nPthatFiles; ++pI) {
    std::string pthatStr = tnamedPthat.at(pI);
    pthatStr.replace(0, pthatStr.find("Pthat") + 5, "");
    pthatFiles[pI] = std::stod(pthatStr);
    TNamed *temp = (TNamed *) inStagFile_p->Get(tnamedPthat.at(pI).c_str());
    std::string xSectionStr = temp->GetTitle();
    xSections[pI] = std::stod(xSectionStr);
  }
  pthatFiles[nPthatFiles] = 9999.;
  xSections[nPthatFiles] = 0.000;

  //do sort
  int pos = 0;
  while (pos < nPthatFiles) {
    bool isGood = true;

    for (int pI = pos + 1; pI < nPthatFiles + 1; ++pI) {
      if (pthatFiles[pos] > pthatFiles[pI]) {
        Double_t tempPthat = pthatFiles[pI];
        Double_t tempX = xSections[pI];

        pthatFiles[pI] = pthatFiles[pos];
        xSections[pI] = xSections[pos];

        pthatFiles[pos] = tempPthat;
        xSections[pos] = tempX;

        isGood = false;
      }
    }

    if (isGood)
      ++pos;
  }

  std::cout << "Sorted xsections: " << std::endl;
  for (Int_t pI = 0; pI < nPthatFiles + 1; ++pI) {
    std::cout << " " << pI << "/" << nPthatFiles + 1 << ": " << pthatFiles[pI] << ", " << xSections[pI] << std::endl;
  }

  inStagFile_p->Close();
  delete inStagFile_p;
  inStagFile_p = NULL;

  Double_t nEvtPerPthatStag[nPthatFiles];
  Double_t nEvtPerPthatStag1200 = 0.;
  Double_t weightsPerPthat[nPthatFiles];

  Double_t nEvtPerPthatFlat[nPthatFiles];

  for (Int_t pI = 0; pI < nPthatFiles; ++pI) {
    nEvtPerPthatStag[pI] = 0.;
    nEvtPerPthatFlat[pI] = 0.;
    weightsPerPthat[pI] = 0.;
  }


  // note: use uniform binning
  // inclusive jet pt spectrum
  TH1D *jetpt_flatPthat_Unweighted_h = new TH1D("jetpt_flatPthat_Unweighted_h",";p_{T} (Flat Unweighted); Counts (Unweighted)", 100, 0, 1000);
  TH1D *jetpt_flatPthat_Weighted_h = new TH1D("jetpt_flatPthat_Weighted_h",";p_{T} (Flat Weighted); Counts (Weighted)", 100, 0, 1000);
  TH1D *jetpt_stagPthat_Unweighted_h = new TH1D("jetpt_stagPthat_Unweighted_h",";p_{T} (Stag Unweighted); Counts (Unweighted)", 100, 0, 1000);
  TH1D *jetpt_stagPthat_Weighted_h = new TH1D("jetpt_stagPthat_Weighted_h",";p_{T} (Stag Weighted); Counts (Weighted)", 100, 0, 1000);


  // dijet xj
  TH1D *dijet_flatPthat_Unweighted_h = new TH1D("dijet_flatPthat_Unweighted_h",
                                                ";X_{J} (Flat Unweighted);Counts (Unweighted)", 10, 0, 1);
  TH1D *dijet_flatPthat_Weighted_h = new TH1D("dijet_flatPthat_Weighted_h",
                                              ";X_{J} (Flat Weighted);Counts (Weighted)", 10, 0, 1);
  TH1D *dijet_stagPthat_Unweighted_h = new TH1D("dijet_stagPthat_Unweighted_h",
                                                ";X_{J} (Stag Unweighted);Counts (Unweighted)", 10, 0, 1);
  TH1D *dijet_stagPthat_Weighted_h = new TH1D("dijet_stagPthat_Weighted_h",
                                              ";X_{J} (Stag Weighted);Counts (Weighted)", 10, 0, 1);
  //dijet dphi
  TH1D *dphi_flatPthat_Unweighted_h = new TH1D("dphi_flatPthat_Unweighted_h",
                                               ";|#Delta#phi| (Flat Unweighted);Counts (Unweighted)", 20, 0, PI);
  TH1D *dphi_flatPthat_Weighted_h = new TH1D("dphi_flatPthat_Weighted_h",
                                             ";|#Delta#phi| (Flat Weighted);Counts (Weighted)", 20, 0, PI);
  TH1D *dphi_stagPthat_Unweighted_h = new TH1D("dphi_stagPthat_Unweighted_h",
                                               ";|#Delta#phi| (Stag Unweighted);Counts (Unweighted)", 20, 0, PI);
  TH1D *dphi_stagPthat_Weighted_h = new TH1D("dphi_stagPthat_Weighted_h",
                                             ";|#Delta#phi| (Stag Weighted);Counts (Weighted)", 20, 0, PI);
  std::vector < TH1 * >tempVect = {
  jetpt_flatPthat_Unweighted_h,jetpt_flatPthat_Weighted_h,
  jetpt_stagPthat_Unweighted_h,jetpt_stagPthat_Weighted_h,
  dijet_flatPthat_Unweighted_h, dijet_flatPthat_Weighted_h,
      dijet_stagPthat_Unweighted_h, dijet_stagPthat_Weighted_h,
      dphi_flatPthat_Unweighted_h, dphi_flatPthat_Weighted_h, dphi_stagPthat_Unweighted_h, dphi_stagPthat_Weighted_h};
  centerTitles(tempVect);
  setSumW2(tempVect);

  Float_t pthat_;
  Float_t weight_;
  Int_t n_jet;
  std::vector < Float_t > *jet_pt = new std::vector < Float_t > ();
  std::vector < Float_t > *jet_phi = new std::vector < Float_t > ();
  std::vector < Float_t > *jet_eta = new std::vector < Float_t > ();

  TFile *inFlatFile_p = new TFile(flatPthatFileName.c_str(), "READ");
  TTree *flatGenTree_p = (TTree *) inFlatFile_p->Get("genTree");

  flatGenTree_p->SetBranchStatus("*", 0);
  flatGenTree_p->SetBranchStatus("pthat", 1);
  flatGenTree_p->SetBranchStatus("weight", 1);
  flatGenTree_p->SetBranchStatus("nGenJt", 1);
  flatGenTree_p->SetBranchStatus("genJtPt", 1);
  flatGenTree_p->SetBranchStatus("genJtPhi", 1);
  flatGenTree_p->SetBranchStatus("genJtEta", 1);

  flatGenTree_p->SetBranchAddress("pthat", &pthat_);
  flatGenTree_p->SetBranchAddress("weight", &weight_);
  flatGenTree_p->SetBranchAddress("nGenJt", &n_jet);
  flatGenTree_p->SetBranchAddress("genJtPt", &jet_pt);
  flatGenTree_p->SetBranchAddress("genJtPhi", &jet_phi);
  flatGenTree_p->SetBranchAddress("genJtEta", &jet_eta);

  const Int_t nEntries = flatGenTree_p->GetEntries();
  std::cout << "Processing flatGenTree, nEntries=" << nEntries << "..." << std::endl;
  for (Int_t entry = 0; entry < nEntries; ++entry) {
    flatGenTree_p->GetEntry(entry);

    if (pthat_ < 15.)
      pthat_ = 15.;

    for (Int_t pI = 0; pI < nPthatFiles; ++pI) {
      if (pthat_ >= pthatFiles[pI] && pthat_ < pthatFiles[pI + 1]) {
        ++nEvtPerPthatFlat[pI];
        break;
      }
    }
  }

  Double_t weightRenorm = flatGenTree_p->GetMaximum("weight");
  std::cout << "Processing flatGenTree, nEntries=" << nEntries << "..." << std::endl;
  for (Int_t entry = 0; entry < nEntries; ++entry) {
    flatGenTree_p->GetEntry(entry);
    int leading = 0;
    int subleading = 0;
    Float_t leading_pt = -9999;
    Float_t subleading_pt = -9999;
    if (n_jet == 0 || n_jet == 1)
      continue;
    for (int i = 0; i < n_jet; i++) {
    jetpt_flatPthat_Unweighted_h->Fill(jet_pt->at(i));
    jetpt_flatPthat_Weighted_h->Fill(jet_pt->at(i), weight_ / weightRenorm);
      if (jet_pt->at(i) > leading_pt) {
        if (leading_pt > subleading_pt) {
          subleading_pt = leading_pt;
          subleading = leading;
        }
        leading_pt = jet_pt->at(i);
        leading = i;
      } else {
        if (jet_pt->at(i) > subleading_pt) {
          subleading_pt = jet_pt->at(i);
          subleading = i;
        } else
          continue;
      }
    }
    Float_t dphi = acos(cos(jet_phi->at(leading) - jet_phi->at(subleading)));
    dphi_flatPthat_Unweighted_h->Fill(dphi);
    dphi_flatPthat_Weighted_h->Fill(dphi, weight_ / weightRenorm);
    if (dphi < PI * 2. / 3. || jet_pt->at(leading) < 120 || jet_pt->at(subleading) < 30
        || std::fabs(jet_eta->at(leading)) > 2 || std::fabs(jet_eta->at(subleading)) > 2)
      continue;
    Float_t xj = jet_pt->at(subleading) / jet_pt->at(leading);
    dijet_flatPthat_Unweighted_h->Fill(xj);
    dijet_flatPthat_Weighted_h->Fill(xj, weight_ / weightRenorm);
  }

  inFlatFile_p->Close();
  delete inFlatFile_p;

  std::cout << "REOPEN FILE?" << std::endl;
  inStagFile_p = new TFile(stagPthatFileName.c_str(), "READ");
  TTree *stagGenTree_p = (TTree *) inStagFile_p->Get("genTree");

  stagGenTree_p->SetBranchStatus("*", 0);
  stagGenTree_p->SetBranchStatus("pthat", 1);
  stagGenTree_p->SetBranchStatus("weight", 1);
  stagGenTree_p->SetBranchStatus("nGenJt", 1);
  stagGenTree_p->SetBranchStatus("genJtPt", 1);
  stagGenTree_p->SetBranchStatus("genJtPhi", 1);
  stagGenTree_p->SetBranchStatus("genJtEta", 1);

  stagGenTree_p->SetBranchAddress("pthat", &pthat_);
  stagGenTree_p->SetBranchAddress("weight", &weight_);
  stagGenTree_p->SetBranchAddress("nGenJt", &n_jet);
  stagGenTree_p->SetBranchAddress("genJtPt", &jet_pt);
  stagGenTree_p->SetBranchAddress("genJtPhi", &jet_phi);
  stagGenTree_p->SetBranchAddress("genJtEta", &jet_eta);

  const Int_t nEntriesStag = stagGenTree_p->GetEntries();
  std::cout << "Processing stagGenTree, nEntries=" << nEntriesStag << "..." << std::endl;
  for (Int_t entry = 0; entry < nEntriesStag; ++entry) {
    stagGenTree_p->GetEntry(entry);

    if (pthat_ < 15.)
      pthat_ = 15.;
    if (pthat_ >= 1200.)
      ++nEvtPerPthatStag1200;

    for (Int_t pI = 0; pI < nPthatFiles; ++pI) {
      if (pthat_ >= pthatFiles[pI] && pthat_ < pthatFiles[pI + 1]) {
        ++nEvtPerPthatStag[pI];
        break;
      }
    }
  }

  std::cout << "Number / pthat: Staggered, Flat (not rescaled)" << std::endl;
  for (Int_t pI = 0; pI < nPthatFiles; ++pI) {
    std::cout << " " << pI << "/" << nPthatFiles << " (" << pthatFiles[pI] << "): "
      << nEvtPerPthatStag[pI] << ", " << nEvtPerPthatFlat[pI] << std::endl;

    weightsPerPthat[pI] = (xSections[pI] - xSections[pI + 1]) / (nEvtPerPthatStag[pI]);
  }

  std::cout << "Number / pthat: Staggered, Flat (rescaled to a 1400000 req.)" << std::endl;
  for (Int_t pI = 0; pI < nPthatFiles; ++pI) {
    std::cout << " " << pI << "/" << nPthatFiles << " (" << pthatFiles[pI] << "): "
      << nEvtPerPthatStag[pI] * 1400000. /
      nEntriesStag << ", " << nEvtPerPthatFlat[pI] * 1400000. / nEntries << std::endl;
  }

  std::cout << "Fraction > 800 that comes from > 1200 in staggered: " << std::endl;
  std::cout << " " << nEvtPerPthatStag1200 << "/" << nEvtPerPthatStag[nPthatFiles -
                                                                      1] << "=" <<
    nEvtPerPthatStag1200 / nEvtPerPthatStag[nPthatFiles - 1] << std::endl;

  for (Int_t pI = 1; pI < nPthatFiles; ++pI) {
    weightsPerPthat[pI] /= weightsPerPthat[0];
  }
  weightsPerPthat[0] = 1.;

  std::cout << "Processing stagGenTree, nEntries=" << nEntriesStag << "..." << std::endl;
  for (Int_t entry = 0; entry < nEntriesStag; ++entry) {
    stagGenTree_p->GetEntry(entry);

    if (pthat_ < 15.)
      pthat_ = 15.;

    Double_t tempWeight_ = 1.;
    for (Int_t pI = 0; pI < nPthatFiles; ++pI) {
      if (pthat_ >= pthatFiles[pI] && pthat_ < pthatFiles[pI + 1]) {
        tempWeight_ = weightsPerPthat[pI];
        break;
      }
    }
    int leading = 0;
    int subleading = 0;
    Float_t leading_pt = -9999;
    Float_t subleading_pt = -9999;
    if (n_jet == 0 || n_jet == 1)
      continue;
    for (int i = 0; i < n_jet; i++) {
    jetpt_stagPthat_Unweighted_h->Fill(jet_pt->at(i));
    jetpt_stagPthat_Weighted_h->Fill(jet_pt->at(i), tempWeight_);
      if (jet_pt->at(i) > leading_pt) {
        if (leading_pt > subleading_pt) {
          subleading_pt = leading_pt;
          subleading = leading;
        }
        leading_pt = jet_pt->at(i);
        leading = i;
      } else {
        if (jet_pt->at(i) > subleading_pt) {
          subleading_pt = jet_pt->at(i);
          subleading = i;
        } else
          continue;
      }
    }
    Float_t dphi = acos(cos(jet_phi->at(leading) - jet_phi->at(subleading)));
    dphi_stagPthat_Unweighted_h->Fill(dphi);
    dphi_stagPthat_Weighted_h->Fill(dphi, tempWeight_);
    if (dphi < PI * 2. / 3. || jet_pt->at(leading) < 120 || jet_pt->at(subleading) < 30
        || std::fabs(jet_eta->at(leading)) > 2 || std::fabs(jet_eta->at(subleading)) > 2)
      continue;
    Float_t xj = jet_pt->at(subleading) / jet_pt->at(leading);
    dijet_stagPthat_Unweighted_h->Fill(xj);
    dijet_stagPthat_Weighted_h->Fill(xj, tempWeight_);
  }

  inStagFile_p->Close();
  delete inStagFile_p;

  outFile_p->cd();
  jetpt_flatPthat_Unweighted_h->Write("", TObject::kOverwrite);
  jetpt_flatPthat_Weighted_h->Write("", TObject::kOverwrite);
  jetpt_stagPthat_Unweighted_h->Write("", TObject::kOverwrite);
  jetpt_stagPthat_Weighted_h->Write("", TObject::kOverwrite);
  dijet_flatPthat_Unweighted_h->Write("", TObject::kOverwrite);
  dijet_flatPthat_Weighted_h->Write("", TObject::kOverwrite);
  dphi_flatPthat_Unweighted_h->Write("", TObject::kOverwrite);
  dphi_flatPthat_Weighted_h->Write("", TObject::kOverwrite);
  dijet_stagPthat_Unweighted_h->Write("", TObject::kOverwrite);
  dijet_stagPthat_Weighted_h->Write("", TObject::kOverwrite);
  dphi_stagPthat_Unweighted_h->Write("", TObject::kOverwrite);
  dphi_stagPthat_Weighted_h->Write("", TObject::kOverwrite);

  //Do Clones
  TH1D *dijet_flatPthat_Unweighted_Scale_h =
    (TH1D *) dijet_flatPthat_Unweighted_h->Clone("dijet_flatPthat_Unweighted_Scale_h");
  TH1D *dijet_flatPthat_Unweighted_ScaleRelErr_h =
    (TH1D *) dijet_flatPthat_Unweighted_h->Clone("dijet_flatPthat_Unweighted_ScaleRelErr_h");
  TH1D *dphi_flatPthat_Unweighted_Scale_h =
    (TH1D *) dphi_flatPthat_Unweighted_h->Clone("dphi_flatPthat_Unweighted_Scale_h");
  TH1D *dphi_flatPthat_Unweighted_ScaleRelErr_h =
    (TH1D *) dphi_flatPthat_Unweighted_h->Clone("dphi_flatPthat_Unweighted_ScaleRelErr_h");

  dijet_flatPthat_Unweighted_ScaleRelErr_h->GetYaxis()->SetTitle("Relative Error, 140k stat");
  dphi_flatPthat_Unweighted_ScaleRelErr_h->GetYaxis()->SetTitle("Relative Error, 140k stat");

  std::cout << "dijet_flatPthat scale" << std::endl;
  Double_t flatScaleFact = 1400000. / nEntries;
  for (Int_t bIX = 0; bIX < dijet_flatPthat_Unweighted_Scale_h->GetNbinsX(); ++bIX) {
    Double_t binVal = dijet_flatPthat_Unweighted_Scale_h->GetBinContent(bIX + 1) * flatScaleFact;
    Double_t binErr = dijet_flatPthat_Unweighted_Scale_h->GetBinError(bIX + 1) * flatScaleFact;
    dijet_flatPthat_Unweighted_Scale_h->SetBinContent(bIX + 1, binVal);
    dijet_flatPthat_Unweighted_Scale_h->SetBinError(bIX + 1, binErr);

    Double_t binRelErr = 0.0;
    Double_t binRelErrErr = 0.0;

    if (binVal > 1) {
      binRelErr = 1. / TMath::Sqrt(binVal);
      binRelErrErr =
        TMath::Max(TMath::Abs(1. / TMath::Sqrt(binVal - binErr) - binRelErr),
                   TMath::Abs(1. / TMath::Sqrt(binVal + binErr) - binRelErr));
    }

    dijet_flatPthat_Unweighted_ScaleRelErr_h->SetBinContent(bIX + 1, binRelErr);
    dijet_flatPthat_Unweighted_ScaleRelErr_h->SetBinError(bIX + 1, binRelErrErr);
  }

  dijet_flatPthat_Unweighted_Scale_h->Write("", TObject::kOverwrite);
  dijet_flatPthat_Unweighted_ScaleRelErr_h->Write("", TObject::kOverwrite);

  std::cout << "dphi_flatPthat scale" << std::endl;
  for (Int_t bIX = 0; bIX < dphi_flatPthat_Unweighted_Scale_h->GetNbinsX(); ++bIX) {
    Double_t binVal = dphi_flatPthat_Unweighted_Scale_h->GetBinContent(bIX + 1) * flatScaleFact;
    Double_t binErr = dphi_flatPthat_Unweighted_Scale_h->GetBinError(bIX + 1) * flatScaleFact;
    dphi_flatPthat_Unweighted_Scale_h->SetBinContent(bIX + 1, binVal);
    dphi_flatPthat_Unweighted_Scale_h->SetBinError(bIX + 1, binErr);

    Double_t binRelErr = 0.0;
    Double_t binRelErrErr = 0.0;

    if (binVal > 1) {
      binRelErr = 1. / TMath::Sqrt(binVal);
      binRelErrErr =
        TMath::Max(TMath::Abs(1. / TMath::Sqrt(binVal - binErr) - binRelErr),
                   TMath::Abs(1. / TMath::Sqrt(binVal + binErr) - binRelErr));
    }

    dphi_flatPthat_Unweighted_ScaleRelErr_h->SetBinContent(bIX + 1, binRelErr);
    dphi_flatPthat_Unweighted_ScaleRelErr_h->SetBinError(bIX + 1, binRelErrErr);
  }

  dphi_flatPthat_Unweighted_Scale_h->Write("", TObject::kOverwrite);
  dphi_flatPthat_Unweighted_ScaleRelErr_h->Write("", TObject::kOverwrite);

  TH1D *dijet_flatPthat_Weighted_Norm_h = (TH1D *) dijet_flatPthat_Weighted_h->Clone("dijet_flatPthat_Weighted_Norm_h");
  dijet_flatPthat_Weighted_Norm_h->Scale(1. / dijet_flatPthat_Weighted_h->Integral());
  dijet_flatPthat_Weighted_Norm_h->Write("", TObject::kOverwrite);
  TH1D *dphi_flatPthat_Weighted_Norm_h = (TH1D *) dphi_flatPthat_Weighted_h->Clone("dphi_flatPthat_Weighted_Norm_h");
  dphi_flatPthat_Weighted_Norm_h->Scale(1. / dphi_flatPthat_Weighted_h->Integral());
  dphi_flatPthat_Weighted_Norm_h->Write("", TObject::kOverwrite);


  std::cout << "dijet_stagPthat scale" << std::endl;
  TH1D *dijet_stagPthat_Weighted_Norm_h = (TH1D *) dijet_stagPthat_Weighted_h->Clone("dijet_stagPthat_Weighted_Norm_h");
  dijet_stagPthat_Weighted_Norm_h->Scale(1. / dijet_stagPthat_Weighted_h->Integral());
  dijet_stagPthat_Weighted_Norm_h->Write("", TObject::kOverwrite);
  TH1D *dphi_stagPthat_Weighted_Norm_h = (TH1D *) dphi_stagPthat_Weighted_h->Clone("dphi_stagPthat_Weighted_Norm_h");
  dphi_stagPthat_Weighted_Norm_h->Scale(1. / dphi_stagPthat_Weighted_h->Integral());
  dphi_stagPthat_Weighted_Norm_h->Write("", TObject::kOverwrite);

  //Do Clones
  TH1D *dijet_stagPthat_Unweighted_Scale_h =
    (TH1D *) dijet_stagPthat_Unweighted_h->Clone("dijet_stagPthat_Unweighted_Scale_h");
  TH1D *dijet_stagPthat_Unweighted_ScaleRelErr_h =
    (TH1D *) dijet_stagPthat_Unweighted_h->Clone("dijet_stagPthat_Unweighted_ScaleRelErr_h");

  dijet_stagPthat_Unweighted_ScaleRelErr_h->GetYaxis()->SetTitle("Relative Error, 1400000 stat");

  Double_t stagScaleFact = 1400000. / nEntriesStag;
  for (Int_t bIX = 0; bIX < dijet_stagPthat_Unweighted_Scale_h->GetNbinsX(); ++bIX) {
    Double_t binVal = dijet_stagPthat_Unweighted_Scale_h->GetBinContent(bIX + 1) * stagScaleFact;
    Double_t binErr = dijet_stagPthat_Unweighted_Scale_h->GetBinError(bIX + 1) * stagScaleFact;
    dijet_stagPthat_Unweighted_Scale_h->SetBinContent(bIX + 1, binVal);
    dijet_stagPthat_Unweighted_Scale_h->SetBinError(bIX + 1, binErr);

    Double_t binRelErr = 0.0;
    Double_t binRelErrErr = 0.0;

    if (binVal > 1) {
      binRelErr = 1. / TMath::Sqrt(binVal);
      binRelErrErr =
        TMath::Max(TMath::Abs(1. / TMath::Sqrt(binVal - binErr) - binRelErr),
                   TMath::Abs(1. / TMath::Sqrt(binVal + binErr) - binRelErr));
    }

    dijet_stagPthat_Unweighted_ScaleRelErr_h->SetBinContent(bIX + 1, binRelErr);
    dijet_stagPthat_Unweighted_ScaleRelErr_h->SetBinError(bIX + 1, binRelErrErr);
  }

  std::cout << "dphi_stagPthat scale" << std::endl;
  dijet_stagPthat_Unweighted_Scale_h->Write("", TObject::kOverwrite);
  dijet_stagPthat_Unweighted_ScaleRelErr_h->Write("", TObject::kOverwrite);

  TH1D *dphi_stagPthat_Unweighted_Scale_h =
    (TH1D *) dphi_stagPthat_Unweighted_h->Clone("dphi_stagPthat_Unweighted_Scale_h");
  TH1D *dphi_stagPthat_Unweighted_ScaleRelErr_h =
    (TH1D *) dphi_stagPthat_Unweighted_h->Clone("dphi_stagPthat_Unweighted_ScaleRelErr_h");

  dphi_stagPthat_Unweighted_ScaleRelErr_h->GetYaxis()->SetTitle("Relative Error, 1400000 stat");

  for (Int_t bIX = 0; bIX < dphi_stagPthat_Unweighted_Scale_h->GetNbinsX(); ++bIX) {
    Double_t binVal = dphi_stagPthat_Unweighted_Scale_h->GetBinContent(bIX + 1) * stagScaleFact;
    Double_t binErr = dphi_stagPthat_Unweighted_Scale_h->GetBinError(bIX + 1) * stagScaleFact;
    dphi_stagPthat_Unweighted_Scale_h->SetBinContent(bIX + 1, binVal);
    dphi_stagPthat_Unweighted_Scale_h->SetBinError(bIX + 1, binErr);

    Double_t binRelErr = 0.0;
    Double_t binRelErrErr = 0.0;

    if (binVal > 1) {
      binRelErr = 1. / TMath::Sqrt(binVal);
      binRelErrErr =
        TMath::Max(TMath::Abs(1. / TMath::Sqrt(binVal - binErr) - binRelErr),
                   TMath::Abs(1. / TMath::Sqrt(binVal + binErr) - binRelErr));
    }

    dphi_stagPthat_Unweighted_ScaleRelErr_h->SetBinContent(bIX + 1, binRelErr);
    dphi_stagPthat_Unweighted_ScaleRelErr_h->SetBinError(bIX + 1, binRelErrErr);
  }

  dphi_stagPthat_Unweighted_Scale_h->Write("", TObject::kOverwrite);
  dphi_stagPthat_Unweighted_ScaleRelErr_h->Write("", TObject::kOverwrite);

  delete dijet_flatPthat_Unweighted_h;
  delete dijet_flatPthat_Unweighted_Scale_h;
  delete dijet_flatPthat_Unweighted_ScaleRelErr_h;
  delete dijet_flatPthat_Weighted_h;
  delete dijet_flatPthat_Weighted_Norm_h;

  delete dijet_stagPthat_Unweighted_h;
  delete dijet_stagPthat_Weighted_h;
  delete dijet_stagPthat_Weighted_Norm_h;
  delete dijet_stagPthat_Unweighted_Scale_h;
  delete dijet_stagPthat_Unweighted_ScaleRelErr_h;

  delete dphi_flatPthat_Unweighted_h;
  delete dphi_flatPthat_Unweighted_Scale_h;
  delete dphi_flatPthat_Unweighted_ScaleRelErr_h;
  delete dphi_flatPthat_Weighted_h;
  delete dphi_flatPthat_Weighted_Norm_h;

  delete dphi_stagPthat_Unweighted_h;
  delete dphi_stagPthat_Weighted_h;
  delete dphi_stagPthat_Weighted_Norm_h;
  delete dphi_stagPthat_Unweighted_Scale_h;
  delete dphi_stagPthat_Unweighted_ScaleRelErr_h;
  
  delete jetpt_flatPthat_Unweighted_h;
  delete jetpt_flatPthat_Weighted_h;
  delete jetpt_stagPthat_Unweighted_h;
  delete jetpt_stagPthat_Weighted_h;


  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: ./bin/statisticalComparison.exe <flatPthatFileName> <stagPthatFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += DijetImbalanceRatio(argv[1], argv[2]);
  return retVal;
}
