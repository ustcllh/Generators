#include "PYTHIA8/PYTHIA.h"

#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>


int pythia8Test(int nEvent)
{
  TFile* outFile_p = new TFile("testWFile.root", "UPDATE");
  TTree* outTree_p = new TTree("genWTree", "genWTree");

  const Int_t genMax = 10000;
  Int_t nGen_;
  Int_t genID_[genMax], genChg_[genMax];
  Float_t genPt_[genMax], genPhi_[genMax], genEta_[genMax];

  Float_t wStartPt_, wStartPhi_, wStartEta_, wStartPx_, wStartPy_;
  Float_t wEndPt_, wEndPhi_, wEndEta_, wEndPx_, wEndPy_;

  Float_t leptGenPt_, leptGenPhi_, leptGenEta_, leptGenPx_, leptGenPy_;
  Float_t neuGenPt_, neuGenPhi_, neuGenEta_, neuGenPx_, neuGenPy_;

  Bool_t isEle_, isMu_;

  outTree_p->Branch("nGen", &nGen_, "nGen/I");
  outTree_p->Branch("genID", genID_, "genID[nGen]/I");
  outTree_p->Branch("genChg", genChg_, "genChg[nGen]/I");
  outTree_p->Branch("genPt", genPt_, "genPt[nGen]/F");
  outTree_p->Branch("genPhi", genPhi_, "genPhi[nGen]/F");
  outTree_p->Branch("genEta", genEta_, "genEta[nGen]/F");

  outTree_p->Branch("wStartPt", &wStartPt_, "wStartPt/F");
  outTree_p->Branch("wStartPhi", &wStartPhi_, "wStartPhi/F");
  outTree_p->Branch("wStartEta", &wStartEta_, "wStartEta/F");
  outTree_p->Branch("wStartPx", &wStartPx_, "wStartPx/F");
  outTree_p->Branch("wStartPy", &wStartPy_, "wStartPy/F");

  outTree_p->Branch("wEndPt", &wEndPt_, "wEndPt/F");
  outTree_p->Branch("wEndPhi", &wEndPhi_, "wEndPhi/F");
  outTree_p->Branch("wEndEta", &wEndEta_, "wEndEta/F");
  outTree_p->Branch("wEndPx", &wEndPx_, "wEndPx/F");
  outTree_p->Branch("wEndPy", &wEndPy_, "wEndPy/F");

  outTree_p->Branch("leptGenPt", &leptGenPt_, "leptGenPt/F");
  outTree_p->Branch("leptGenPhi", &leptGenPhi_, "leptGenPhi/F");
  outTree_p->Branch("leptGenEta", &leptGenEta_, "leptGenEta/F");
  outTree_p->Branch("leptGenPx", &leptGenPx_, "leptGenPx/F");
  outTree_p->Branch("leptGenPy", &leptGenPy_, "leptGenPy/F");

  outTree_p->Branch("neuGenPt", &neuGenPt_, "neuGenPt/F");
  outTree_p->Branch("neuGenPhi", &neuGenPhi_, "neuGenPhi/F");
  outTree_p->Branch("neuGenEta", &neuGenEta_, "neuGenEta/F");
  outTree_p->Branch("neuGenPx", &neuGenPx_, "neuGenPx/F");
  outTree_p->Branch("neuGenPy", &neuGenPy_, "neuGenPy/F");

  outTree_p->Branch("isEle", &isEle_, "isEle/O");
  outTree_p->Branch("isMu", &isMu_, "isMu/O");

  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("WeakBosonAndParton:qqbar2Wg = on");
  pythia.readString("WeakBosonAndParton:qg2Wq = on");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = 11 12 13 14");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("PhaseSpace:pTHatMin = 80");


  pythia.init();

  for(int iter = 0; iter < nEvent*10; iter++){
    if(iter%1000 == 0) std::cout << "Event, Fills: " << iter << ", " << outTree_p->GetEntries() << std::endl;

    nGen_ = 0;
    wStartPt_ = -999;
    wStartPhi_ = -999;
    wStartEta_ = -999;
    wStartPx_ = -999;
    wStartPy_ = -999;

    wEndPt_ = -999;
    wEndPhi_ = -999;
    wEndEta_ = -999;
    wEndPx_ = -999;
    wEndPy_ = -999;

    leptGenPt_ = -999;
    leptGenPhi_ = -999;
    leptGenEta_ = -999;
    leptGenPx_ = -999;
    leptGenPy_ = -999;

    neuGenPt_ = -999;
    neuGenPhi_ = -999;
    neuGenEta_ = -999;
    neuGenPx_ = -999;
    neuGenPy_ = -999;

    isEle_ = false;
    isMu_ = false;

    pythia.next();

    for(int entry = 0; entry < pythia.event.size(); entry++){
      if(std::abs(pythia.event[entry].id()) == 24){
	if(wStartPt_ < 0){
	  wStartPt_ = pythia.event[entry].pT();
	  wStartPhi_ = pythia.event[entry].phi();
	  wStartEta_ = pythia.event[entry].eta();
	  wStartPx_ = pythia.event[entry].px();
	  wStartPy_ = pythia.event[entry].py();
	}
	else{
	  wEndPt_ = pythia.event[entry].pT();
	  wEndPhi_ = pythia.event[entry].phi();
	  wEndEta_ = pythia.event[entry].eta();
	  wEndPx_ = pythia.event[entry].px();
	  wEndPy_ = pythia.event[entry].py();
	}
      }
      else if(std::abs(pythia.event[entry].id()) == 11 || std::abs(pythia.event[entry].id()) == 13){
	if(std::abs(pythia.event[pythia.event[entry].mother1()].id()) == 24 && std::abs(pythia.event[entry].mother2()) == 0){
	  leptGenPt_ = pythia.event[entry].pT();
          leptGenPhi_ = pythia.event[entry].phi();
          leptGenEta_ = pythia.event[entry].eta();
	  leptGenPx_ = pythia.event[entry].px();
	  leptGenPy_ = pythia.event[entry].py();

	  if(TMath::Abs(leptGenEta_) > 1.44 && leptGenPt_ > 0) break;

	  if(std::abs(pythia.event[entry].id()) == 11) isEle_ = true;
	  else isMu_ = true;
	}
      }
      else if(std::abs(pythia.event[entry].id()) == 12 || std::abs(pythia.event[entry].id()) == 14){
	if(std::abs(pythia.event[pythia.event[entry].mother1()].id()) == 24 && std::abs(pythia.event[entry].mother2()) == 0){
	  neuGenPt_ = pythia.event[entry].pT();
          neuGenPhi_ = pythia.event[entry].phi();
          neuGenEta_ = pythia.event[entry].eta();
	  neuGenPx_ = pythia.event[entry].px();
	  neuGenPy_ = pythia.event[entry].py();

	  if(TMath::Abs(neuGenEta_) > 2.4 && neuGenPt_ > 0) break;
	}
      }

      if(pythia.event[entry].isFinal()){
	genID_[nGen_] = pythia.event[entry].id();
	genChg_[nGen_] = pythia.event[entry].charge();
	genPt_[nGen_] = pythia.event[entry].pT();
	genPhi_[nGen_] = pythia.event[entry].phi();
	genEta_[nGen_] = pythia.event[entry].eta();
	nGen_++;
      }
    }

    if(TMath::Abs(leptGenEta_) > 1.44 && leptGenPt_ > 0) continue;
    if(TMath::Abs(neuGenEta_) > 2.4 && neuGenPt_ > 0) continue;

    outTree_p->Fill();

    if(outTree_p->GetEntries() == nEvent) break;
  }

  outFile_p->cd();
  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: pythia8Test <nEvents>" << std::endl;

    std::cout << "argc: " << argc << std::endl;

    for(int iter = 0; iter < argc; iter++){
      std::cout << "  argv[" << iter << "]: " << argv[iter] << std::endl;
    }

    return 1;
  }

  int rStatus = -1;
  rStatus = pythia8Test(std::atoi(argv[1]));
  return rStatus;
}
