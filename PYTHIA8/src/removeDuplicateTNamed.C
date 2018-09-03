//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"

//Non-Local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"

int removeDuplicateTNamed(const std::string inFileName, std::string outFileName = "")
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is invaled. return 1" << std::endl;  
    return 1;
  }

  if(outFileName.size() == 0){
    outFileName = inFileName.substr(0, inFileName.find(".root"));
    outFileName = outFileName + "_RemovedDupTNamed.root";
  }


  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> classList;
  std::vector<std::string> contentsList = returnRootFileContentsList(inFile_p, "", "", 1, &(classList));

  std::cout << "Printing contents: " << std::endl;
  for(unsigned int i = 0; i < contentsList.size(); ++i){
    if(isStrSame(classList.at(i), "TTree")){
      
    }
  }

  
  
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/removeDuplicateTNamed.exe <inFileName> <outFileName-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += removeDuplicateTNamed(argv[1]);
  else if(argc == 3) retVal += removeDuplicateTNamed(argv[1], argv[2]);
  return retVal;
}
