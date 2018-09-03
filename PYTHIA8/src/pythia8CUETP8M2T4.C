//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDatime.h"
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TDirectoryFile.h"

//PYTHIA8 dependencies
#include "Pythia8/Pythia.h"

//FASTJET dependencies
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/stringUtil.h"

//TUNE LINKS:
//Common Block: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CommonSettings_cfi.py
//CUETP8M1: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CUEP8M1Settings_cfi.py
//CP5: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/MCTunes2017/PythiaCP5Settings_cfi.py
//CUETP8M2T4: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CUEP8M2T4Settings_cfi.py

int pythia8CUETP8M2T4(std::string outFileName, const std::string tuneStr, const int nEvt = 10000, const bool doFlatPthat = false, const double flatPthatScale = 4.5, const double pthatMin = 15., const double pthatMax = 999999.)
{
  const Int_t nTunes = 4;
  const std::string tunes[nTunes] = {"Vanilla", "CUETP8M1", "CP5", "CUETP8M2T4"};

  int tunePos = -1;
  
  for(Int_t tI = 0; tI < nTunes; ++tI){
    if(isStrSame(tuneStr, tunes[tI])){
      tunePos = tI;
      break;
    }
  }

  if(tunePos < 0){
    std::cout << "Given tuneStr, \'" << tuneStr << "\' is invalid. Please pick one of the following: " << std::endl;
    for(Int_t tI = 0; tI < nTunes; ++tI){
      std::cout << " " << tunes[tI] << ",";
    }
    std::cout << std::endl;
    std::cout << "Return 1" << std::endl;
    return 1;
  }
    

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1,"");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string descStr = tuneStr + "_NEvt" + std::to_string(nEvt) + "_DoFlatPthat";
  if(doFlatPthat) descStr = descStr + "True_FlatPower" + prettyString(flatPthatScale, 2, true);
  else descStr = descStr + "False";
  
  descStr = descStr + "_PthatMin" + prettyString(pthatMin, 1, true) + "_PthatMax" + prettyString(pthatMax, 1, true);
  
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + descStr + "_" + dateStr + ".root";
  
  std::cout << "Creating " << outFileName << "..." << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* genTree_p = new TTree("genTree", "");
  TDirectoryFile* paramsDir_p = (TDirectoryFile*)outFile_p->mkdir("paramsDir");
  
  Float_t pthat_;
  Float_t weight_;

  //define jet parameters
  const float minJtPtCut = 15.;
  const float maxJtAbsEtaCut = 5.1;

  const fastjet::JetAlgorithm jtAlgo = fastjet::antikt_algorithm;
  const double jtRVal = 0.4;
  const fastjet::RecombinationScheme jtRecombScheme = fastjet::E_scheme;
  const fastjet::JetDefinition jtDef(jtAlgo, jtRVal, jtRecombScheme);

  Int_t nGenJt_;
  std::vector<float>* genJtPt_p = new std::vector<float>;
  std::vector<float>* genJtPhi_p = new std::vector<float>;
  std::vector<float>* genJtEta_p = new std::vector<float>;
  
  genTree_p->Branch("pthat", &pthat_, "pthat/F");
  genTree_p->Branch("weight", &weight_, "weight/F");  
  genTree_p->Branch("nGenJt", &nGenJt_, "nGenJt/I");
  genTree_p->Branch("genJtPt", &genJtPt_p);
  genTree_p->Branch("genJtPhi", &genJtPhi_p);
  genTree_p->Branch("genJtEta", &genJtEta_p);  
  
  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");

  //Vanilla implement, shared

  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Main:timesAllowErrors = 10000");
  pythia.readString("Check:epTolErr = 0.01");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");
  pythia.readString("SLHA:keepSM = on");
  pythia.readString("SLHA:minMassSM = 1000.");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  if(tunePos == 1){//CUETP8M1
    pythia.readString("Tune:pp 14");
    pythia.readString("Tune:ee 7");
    pythia.readString("MultipartonInteractions:pT0Ref=2.4024");
    pythia.readString("MultipartonInteractions:ecmPow=0.25208");
    pythia.readString("MultipartonInteractions:expPow=1.6");
  }
  else if(tunePos == 2){//CP5
    pythia.readString("Tune:pp 14");
    pythia.readString("Tune:ee 7");
    pythia.readString("MultipartonInteractions:ecmPow=0.03344");
    pythia.readString("PDF:pSet=20");
    pythia.readString("MultipartonInteractions:bProfile=2");
    pythia.readString("MultipartonInteractions:pT0Ref=1.41");
    pythia.readString("MultipartonInteractions:coreRadius=0.7634");
    pythia.readString("MultipartonInteractions:coreFraction=0.63");
    pythia.readString("ColourReconnection:range=5.176");
    pythia.readString("SigmaTotal:zeroAXB=off");
    pythia.readString("SpaceShower:alphaSorder=2");
    pythia.readString("SpaceShower:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSorder=2");
    pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
    pythia.readString("MultipartonInteractions:alphaSorder=2");
    pythia.readString("TimeShower:alphaSorder=2");
    pythia.readString("TimeShower:alphaSvalue=0.118");
  }
  else if(tunePos == 3){//CUETP8M2T4
    pythia.readString("Tune:pp 14");
    pythia.readString("Tune:ee 7");
    pythia.readString("MultipartonInteractions:ecmPow=0.25208");
    pythia.readString("SpaceShower:alphaSvalue=0.1108");
    pythia.readString("PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130");
    pythia.readString("MultipartonInteractions:pT0Ref=2.20e+00");
    pythia.readString("MultipartonInteractions:expPow=1.60e+00");
    pythia.readString("ColourReconnection:range=6.59e+00");
  }
  
  if(doFlatPthat){
    pythia.readString("PhaseSpace:bias2Selection = on");
    pythia.readString("PhaseSpace:bias2SelectionPow = " + std::to_string(flatPthatScale));
    pythia.readString("PhaseSpace:bias2SelectionRef = " + std::to_string(pthatMin));
  }

  pythia.readString("PhaseSpace:pTHatMin = "  + std::to_string(pthatMin));
  pythia.readString("PhaseSpace:pTHatMax = "  + std::to_string(pthatMax));

  pythia.init();

  int totEvt = 0;
  
  while(totEvt < nEvt){
    if(!pythia.next()) continue;
    
    nGenJt_=0;
    genJtPt_p->clear();
    genJtPhi_p->clear();
    genJtEta_p->clear();
    
    pthat_ = pythia.info.pTHat();
    if(doFlatPthat) weight_ = pythia.info.weight();
    else weight_ = 1.;

    //weight_ = 1.;

    std::vector<fastjet::PseudoJet> particles;
    
    for (int i = 0; i < pythia.event.size(); ++i){
      if(!pythia.event[i].isFinal()) continue;
      if(pythia.event[i].pT() < 0.1) continue; // assuming gev
      if(TMath::Abs(pythia.event[i].eta()) > maxJtAbsEtaCut) continue; // assuming rough detector geometry
      //continuing on neutrinos;
      if(TMath::Abs(pythia.event[i].id()) == 12) continue;
      if(TMath::Abs(pythia.event[i].id()) == 14) continue;
      if(TMath::Abs(pythia.event[i].id()) == 16) continue;

      particles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
    }

    fastjet::ClusterSequence* cs = new fastjet::ClusterSequence(particles, jtDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs->inclusive_jets(minJtPtCut));

    for(unsigned int jI = 0; jI < jets.size(); ++jI){
      if(jets.at(jI).eta() > maxJtAbsEtaCut) continue;

      genJtPt_p->push_back(jets.at(jI).pt());
      genJtPhi_p->push_back(jets.at(jI).phi_std());
      genJtEta_p->push_back(jets.at(jI).eta());
      ++nGenJt_;
    }

    particles.clear();
    jets.clear();
    delete cs;
    
    genTree_p->Fill();
    particles.clear();
    
    ++totEvt;
  }

  double sigmapb = pythia.info.sigmaGen() * 1.0E9;
  std::cout << "== Cross section for this run = " <<  sigmapb << " pb" << std::endl;
  
  outFile_p->cd();

  genTree_p->Write("", TObject::kOverwrite);

  paramsDir_p->cd();

  TNamed tuneName("tune", tuneStr.c_str());
  TNamed nEvtName("nEvt", std::to_string(nEvt).c_str());
  TNamed doFlatPthatName("doFlatPthat", std::to_string(doFlatPthat).c_str());
  TNamed flatPthatScaleName("flatPthatScale", prettyString(flatPthatScale, 2, false).c_str());
  TNamed pthatMinName("pthatMin", prettyString(pthatMin, 1, false).c_str());
  TNamed pthatMaxName("pthatMax", prettyString(pthatMax, 1, false).c_str());
  std::string crossSectionString = "crossSection_Pthat" + std::to_string((int)pthatMin);
  TNamed crossSectionName(crossSectionString.c_str(), prettyString(sigmapb, 3, false).c_str());
  
  tuneName.Write("", TObject::kOverwrite);
  nEvtName.Write("", TObject::kOverwrite);
  doFlatPthatName.Write("", TObject::kOverwrite);
  flatPthatScaleName.Write("", TObject::kOverwrite);
  pthatMinName.Write("", TObject::kOverwrite);
  pthatMaxName.Write("", TObject::kOverwrite);
  crossSectionName.Write("", TObject::kOverwrite);
  
  delete genTree_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 3 || argc > 8){
    std::cout << "Usage: ./bin/pythia8CUETP8M2T4.exe <outFileName> <tuneStr> <nEvt-opt> <doFlatPthat-opt> <flatPthatScale-opt> <pthatMin-opt> <pthatMax-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += pythia8CUETP8M2T4(argv[1], argv[2]);
  else if(argc == 4) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]));
  else if(argc == 5) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]));
  else if(argc == 6) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]));
  else if(argc == 7) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]));
  else if(argc == 8) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]));
  return retVal;
}
