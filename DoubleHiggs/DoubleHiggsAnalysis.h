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
#include "TTree.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootProgressBar.h"


using namespace std;


class DoubleHiggsAnalysis
{
public:
  DoubleHiggsAnalysis(const char *inputFile, const char *outputFile);
  ~DoubleHiggsAnalysis();

  TFile* outFile_;
  TTree* tree_passedEvents;
  TChain* chain_;
  ExRootTreeReader *treeReader_;

  Photon *pho1, *pho2;
  Jet *jet1, *jet2;	  

  float xSec_;
  float totalGenEvents_;
  float eventWeight_;

  //variables for the output tree
  int njets_t, nbjets_t,nleptons_t,nelectrons_t,nmuons_t,electron_index_t,muon_index_t;
  float ptPhot1_t, etaPhot1_t, phiPhot1_t;
  float ptPhot2_t, etaPhot2_t, phiPhot2_t;
  float mgg_t, ptJet_t[10], etaJet_t[10], phiJet_t[10];
  float minDeltaRElePho_t,deltaRGammaGamma_t,minDeltaRGammaB_t, deltaRBB_t;
  float ptEle_t,ptMuon_t;
  float ptGG_t,ptBB_t,ptBBGG_t;
  float mbb_t,mggbb_t;
  float eventWeight_t;
  std::pair<int,int> bjet_indexes;
  std::pair<float,float> bjet_pt;


  bool LeptonSelection(TClonesArray *branchElectron,TClonesArray *branchMuon,Photon *pho1, Photon *pho2);
  bool JetSelection(TClonesArray *branchJet, bool looseBtagWP,Photon *pho1, Photon *pho2);
  bool PhotonSelection(TClonesArray *branchPhoton);

  void createOutputTree();
  void setGenEvents(float genevents);
  void setEventWeight();
  void setXsec(float xsec);
  void setOutFile(const char *outputFile);
  void Analyze();

};
