#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1F.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootProgressBar.h"


using namespace std;


class DoubleHiggsAnalysis
{
public:
  DoubleHiggsAnalysis(const char *inputFile);
  ~DoubleHiggsAnalysis();

  TFile* outFile_;
  TChain* chain_;
  ExRootTreeReader *treeReader_;

  bool JetSelection(TClonesArray *branchJet, bool tightBtagWP);
  bool PhotonSelection(TClonesArray *branchPhoton);

  void setOutFile(const char *outputFile);
  void Analyze();

};
