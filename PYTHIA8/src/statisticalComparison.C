//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TH1D.h"

//Non-Local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"

int statisticalComparison(const std::string flatPthatFileName, const std::string stagPthatFileName)
{
  bool goodFlatFile = checkFile(flatPthatFileName);
  bool goodStagFile = checkFile(stagPthatFileName);

  if(!goodFlatFile || !goodStagFile){
    if(!goodFlatFile) std::cout << "Input flat pTHat file \'" << flatPthatFileName << "\' is invalid. return 1" << std::endl;
    if(!goodStagFile) std::cout << "Input stag pTHat file \'" << stagPthatFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  checkMakeDir("output");

  const std::string outFileName = "output/statComp_" + dateStr + ".root";
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  const Int_t nPthatFiles = 14;
  const Double_t pthatFiles[nPthatFiles+1] = {15., 30., 50., 80., 120., 170., 220., 280., 370., 460., 540., 630., 720., 800., 9999.};
  //Cross sections in pb extracted from the 100k processing of each pthat, output in the log, grabbed with grep
  const Double_t xSections[nPthatFiles+1] = {524329708.942, 34679381.584, 4029107.197, 491741.431, 71146.579, 12224.793, 3058.662, 774.598, 140.047, 32.435, 10.110, 2.970, 0.930, 0.347, 0.0000000};  
  Double_t nEvtPerPthatStag[nPthatFiles];
  Double_t nEvtPerPthatStag1200 = 0.;
  Double_t weightsPerPthat[nPthatFiles];

  Double_t nEvtPerPthatFlat[nPthatFiles];

  for(Int_t pI = 0; pI < nPthatFiles; ++pI){
    nEvtPerPthatStag[pI] = 0.;
    nEvtPerPthatFlat[pI] = 0.;
    weightsPerPthat[pI] = 0.;
  }
  
  const Int_t nPthatBins = 200;
  Double_t pthatBins[nPthatBins+1];
  const Double_t pthatLow = 15.;
  const Double_t pthatHi = 1415.;
  getLinBins(pthatLow, pthatHi, nPthatBins, pthatBins);
  
  TH1D* flatPthat_Unweighted_h = new TH1D("flatPthat_Unweighted_h", ";Flat p_{T} Hat (Unweighted);Counts (Unweighted)", nPthatBins, pthatBins);
  TH1D* flatPthat_Weighted_h = new TH1D("flatPthat_Weighted_h", ";Flat p_{T} Hat (Weighted);Counts (Weighted)", nPthatBins, pthatBins);

  TH1D* stagPthat_Unweighted_h = new TH1D("stagPthat_Unweighted_h", ";Stag p_{T} Hat (Unweighted);Counts (Unweighted)", nPthatBins, pthatBins);
  TH1D* stagPthat_Weighted_h = new TH1D("stagPthat_Weighted_h", ";Stag p_{T} Hat (Weighted);Counts (Weighted)", nPthatBins, pthatBins);

  std::vector<TH1*> tempVect = {flatPthat_Unweighted_h, flatPthat_Weighted_h, stagPthat_Unweighted_h, stagPthat_Weighted_h};
  centerTitles(tempVect);
  setSumW2(tempVect);
  
  Float_t pthat_;
  Float_t weight_;
  
  TFile* inFlatFile_p = new TFile(flatPthatFileName.c_str(), "READ");
  TTree* flatGenTree_p = (TTree*)inFlatFile_p->Get("genTree");
  flatGenTree_p->SetBranchAddress("pthat", &pthat_);
  flatGenTree_p->SetBranchAddress("weight", &weight_);

  Double_t weightRenorm = flatGenTree_p->GetMaximum("weight");
  
  const Int_t nEntries = flatGenTree_p->GetEntries();

  std::cout << "Processing flatGenTree, nEntries=" << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    flatGenTree_p->GetEntry(entry);

    if(pthat_ < 15.) pthat_ = 15.;

    for(Int_t pI = 0; pI < nPthatFiles; ++pI){
      if(pthat_ >= pthatFiles[pI] && pthat_ < pthatFiles[pI+1]){
	++nEvtPerPthatFlat[pI];
	break;
      }
    }    

    flatPthat_Unweighted_h->Fill(pthat_);
    flatPthat_Weighted_h->Fill(pthat_, weight_/weightRenorm);
  }
  
  inFlatFile_p->Close();
  delete inFlatFile_p;

  TFile* inStagFile_p = new TFile(stagPthatFileName.c_str(), "READ");
  TTree* stagGenTree_p = (TTree*)inStagFile_p->Get("genTree");
  stagGenTree_p->SetBranchAddress("pthat", &pthat_);

  const Int_t nEntriesStag = stagGenTree_p->GetEntries();

  std::cout << "Processing stagGenTree, nEntries=" << nEntriesStag << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntriesStag; ++entry){
    stagGenTree_p->GetEntry(entry);

    if(pthat_ < 15.) pthat_ = 15.;
    if(pthat_ >= 1200.) ++nEvtPerPthatStag1200;
    
    for(Int_t pI = 0; pI < nPthatFiles; ++pI){
      if(pthat_ >= pthatFiles[pI] && pthat_ < pthatFiles[pI+1]){
	++nEvtPerPthatStag[pI];
	break;
      }
    }    
  }

  std::cout << "Number / pthat: Staggered, Flat (not rescaled)" << std::endl;
  for(Int_t pI = 0; pI < nPthatFiles; ++pI){
    std::cout << " " << pI << "/" << nPthatFiles << " (" << pthatFiles[pI] << "): " << nEvtPerPthatStag[pI] << ", " << nEvtPerPthatFlat[pI] << std::endl;
    
    weightsPerPthat[pI] = (xSections[pI] - xSections[pI+1])/(nEvtPerPthatStag[pI]);
  }

  std::cout << "Number / pthat: Staggered, Flat (rescaled to a 70 million req.)" << std::endl;
  for(Int_t pI = 0; pI < nPthatFiles; ++pI){
    std::cout << " " << pI << "/" << nPthatFiles << " (" << pthatFiles[pI] << "): " << nEvtPerPthatStag[pI]*70000000./nEntriesStag << ", " << nEvtPerPthatFlat[pI]*70000000./nEntries << std::endl;
  }

  std::cout << "Fraction > 800 that comes from > 1200 in staggered: " << std::endl;
  std::cout << " " << nEvtPerPthatStag1200 << "/" << nEvtPerPthatStag[nPthatFiles-1] << "=" << nEvtPerPthatStag1200/nEvtPerPthatStag[nPthatFiles-1] << std::endl;
  
  for(Int_t pI = 1; pI < nPthatFiles; ++pI){
    weightsPerPthat[pI] /= weightsPerPthat[0];
  }
  weightsPerPthat[0] = 1.;

  std::cout << "Processing stagGenTree, nEntries=" << nEntriesStag << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntriesStag; ++entry){
    stagGenTree_p->GetEntry(entry);

    if(pthat_ < 15.) pthat_ = 15.;

    Double_t tempWeight_ = 1.;
    for(Int_t pI = 0; pI < nPthatFiles; ++pI){
      if(pthat_ >= pthatFiles[pI] && pthat_ < pthatFiles[pI+1]){
	tempWeight_ = weightsPerPthat[pI];
	break;
      }
    }

    stagPthat_Unweighted_h->Fill(pthat_);
    stagPthat_Weighted_h->Fill(pthat_, tempWeight_);
  }


  inStagFile_p->Close();
  delete inStagFile_p;

  outFile_p->cd();

  flatPthat_Unweighted_h->Write("", TObject::kOverwrite);
  flatPthat_Weighted_h->Write("", TObject::kOverwrite);

  //Do Clones
  TH1D* flatPthat_Unweighted_Scale_h = (TH1D*)flatPthat_Unweighted_h->Clone("flatPthat_Unweighted_Scale_h");
  TH1D* flatPthat_Unweighted_ScaleRelErr_h = (TH1D*)flatPthat_Unweighted_h->Clone("flatPthat_Unweighted_ScaleRelErr_h");

  flatPthat_Unweighted_ScaleRelErr_h->GetYaxis()->SetTitle("Relative Error, 7E10stat");
  
  Double_t flatScaleFact = 70000000./nEntries;
  for(Int_t bIX = 0; bIX < flatPthat_Unweighted_Scale_h->GetNbinsX(); ++bIX){
    Double_t binVal = flatPthat_Unweighted_Scale_h->GetBinContent(bIX+1)*flatScaleFact;
    Double_t binErr = flatPthat_Unweighted_Scale_h->GetBinError(bIX+1)*flatScaleFact;
    flatPthat_Unweighted_Scale_h->SetBinContent(bIX+1, binVal);
    flatPthat_Unweighted_Scale_h->SetBinError(bIX+1, binErr);

    Double_t binRelErr = 0.0; 
    Double_t binRelErrErr = 0.0;

    if(binVal > 0){
      binRelErr = 1./TMath::Sqrt(binVal);
      binRelErrErr = TMath::Max(TMath::Abs(1./TMath::Sqrt(binVal - binErr) - binRelErr), TMath::Abs(1./TMath::Sqrt(binVal + binErr) - binRelErr));
    }
    
    flatPthat_Unweighted_ScaleRelErr_h->SetBinContent(bIX+1, binRelErr);
    flatPthat_Unweighted_ScaleRelErr_h->SetBinError(bIX+1, binRelErrErr);
  }  

  flatPthat_Unweighted_Scale_h->Write("", TObject::kOverwrite);
  flatPthat_Unweighted_ScaleRelErr_h->Write("", TObject::kOverwrite);
  
  TH1D* flatPthat_Weighted_Norm_h = (TH1D*)flatPthat_Weighted_h->Clone("flatPthat_Weighted_Norm_h");
  flatPthat_Weighted_Norm_h->Scale(1./flatPthat_Weighted_h->Integral());
  flatPthat_Weighted_Norm_h->Write("", TObject::kOverwrite);

  
  stagPthat_Unweighted_h->Write("", TObject::kOverwrite);
  stagPthat_Weighted_h->Write("", TObject::kOverwrite);
  TH1D* stagPthat_Weighted_Norm_h = (TH1D*)stagPthat_Weighted_h->Clone("stagPthat_Weighted_Norm_h");
  stagPthat_Weighted_Norm_h->Scale(1./stagPthat_Weighted_h->Integral());
  stagPthat_Weighted_Norm_h->Write("", TObject::kOverwrite);

  //Do Clones
  TH1D* stagPthat_Unweighted_Scale_h = (TH1D*)stagPthat_Unweighted_h->Clone("stagPthat_Unweighted_Scale_h");
  TH1D* stagPthat_Unweighted_ScaleRelErr_h = (TH1D*)stagPthat_Unweighted_h->Clone("stagPthat_Unweighted_ScaleRelErr_h");

  stagPthat_Unweighted_ScaleRelErr_h->GetYaxis()->SetTitle("Relative Error, 7E10stat");

  Double_t stagScaleFact = 70000000./nEntriesStag;
  for(Int_t bIX = 0; bIX < stagPthat_Unweighted_Scale_h->GetNbinsX(); ++bIX){
    Double_t binVal = stagPthat_Unweighted_Scale_h->GetBinContent(bIX+1)*stagScaleFact;
    Double_t binErr = stagPthat_Unweighted_Scale_h->GetBinError(bIX+1)*stagScaleFact;
    stagPthat_Unweighted_Scale_h->SetBinContent(bIX+1, binVal);
    stagPthat_Unweighted_Scale_h->SetBinError(bIX+1, binErr);

    Double_t binRelErr = 0.0; 
    Double_t binRelErrErr = 0.0;

    if(binVal > 0){
      binRelErr = 1./TMath::Sqrt(binVal);
      binRelErrErr = TMath::Max(TMath::Abs(1./TMath::Sqrt(binVal - binErr) - binRelErr), TMath::Abs(1./TMath::Sqrt(binVal + binErr) - binRelErr));
    }

    stagPthat_Unweighted_ScaleRelErr_h->SetBinContent(bIX+1, binRelErr);
    stagPthat_Unweighted_ScaleRelErr_h->SetBinError(bIX+1, binRelErrErr);
  }  

  stagPthat_Unweighted_Scale_h->Write("", TObject::kOverwrite);
  stagPthat_Unweighted_ScaleRelErr_h->Write("", TObject::kOverwrite);

  delete flatPthat_Unweighted_h;
  delete flatPthat_Unweighted_Scale_h;
  delete flatPthat_Unweighted_ScaleRelErr_h;
  delete flatPthat_Weighted_h;
  delete flatPthat_Weighted_Norm_h;

  delete stagPthat_Unweighted_h;
  delete stagPthat_Weighted_h;
  delete stagPthat_Weighted_Norm_h;
  delete stagPthat_Unweighted_Scale_h;
  delete stagPthat_Unweighted_ScaleRelErr_h;

  outFile_p->Close();
  delete outFile_p;
  
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/statisticalComparison.exe <flatPthatFileName> <stagPthatFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += statisticalComparison(argv[1], argv[2]);
  return retVal;
}
