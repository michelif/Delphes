#include "DoubleHiggsAnalysis.h"

void DoubleHiggsAnalysis::setOutFile(const char *outputFile){

  outFile_=TFile::Open(outputFile,"recreate");

}


void DoubleHiggsAnalysis::Analyze(){

  Long64_t numberOfEntries = treeReader_->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader_->UseBranch("Jet");
  TClonesArray *branchPhoton = treeReader_->UseBranch("Photon");

  TH1F *histPhoMass = new TH1F("mass", "M_{inv}(#gamma#gamma)", 100, 100.0, 180.0);

  setOutFile("outprova.root");
  cout<<"starting loop"<<endl;
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
    // Load selected branches with data from specified event
      treeReader_->ReadEntry(entry);
      bool passedJetSelection=false;
      bool passedPhotonSelection=false;
      
      bool tightBtagWP=false;
      
      passedJetSelection=JetSelection(branchJet, tightBtagWP);
      cout<<"passed jet"<<endl;
      passedPhotonSelection=PhotonSelection(branchPhoton);
      cout<<"passed pho"<<endl;      
      
      //filling histograms
      if(passedJetSelection && passedPhotonSelection){
	Photon *pho1, *pho2;
	
	// Take first two photons
	pho1 = (Photon *) branchPhoton->At(0);
	pho2 = (Photon *) branchPhoton->At(1);
  
      histPhoMass->Fill(((pho1->P4()) + (pho2->P4())).M());
      }
    }
  
  // Show resulting histograms
  //  histJetPT->Write();
  //  histPhoMass->Write();
  //  histPhoLeadPT->Write();
  //  histPhoSubleadPT->Write();
  histPhoMass->Write();
  outFile_->Write();
  outFile_->Close();

}

DoubleHiggsAnalysis::DoubleHiggsAnalysis(const char *inputFile)
{


  // Create chain of root trees
  chain_=new TChain("Delphes");
  chain_->Add(inputFile);
  
  // Create object of class ExRootTreeReader
  treeReader_ = new ExRootTreeReader(chain_);
  

}


//-------------------------------------------
bool DoubleHiggsAnalysis::JetSelection(TClonesArray *branchJet, bool tightBtagWP){

  bool passed=false;

  // If event contains at least 2 jet
  if(branchJet->GetEntries() < 2) return passed;
    
  std::pair<int,int> bjet_indexes;
  std::pair<float,float> bjet_pt;
  
  bjet_indexes.first=-1;
  bjet_indexes.second=-1;
  bjet_pt.first=-1.;
  bjet_pt.second=-1.;
  
  
  for(int ijet=0;ijet<branchJet->GetEntries();ijet++){
    Jet *jet = (Jet*) branchJet->At(ijet);
    if(TMath::Abs(jet->Eta>2.4)) continue;

    bool isBtaggedNormal = (jet->BTag & (1 << 0));
    bool isBtaggedLoose= (jet->BTag & (1 << 1));

    bool isBtagged=isBtaggedLoose;
    if(tightBtagWP)isBtagged=isBtaggedNormal;
    
    if(!isBtagged) continue;

    //bjets
    if(isBtagged){
      if (jet->PT>bjet_pt.first){
	bjet_indexes.second=bjet_indexes.first;
	bjet_pt.second=bjet_pt.first;
	bjet_indexes.first=ijet;
	bjet_pt.first=jet->PT;
      }else if (jet->PT>bjet_pt.second){
	bjet_indexes.second=ijet;
	bjet_pt.second=jet->PT;
      }
    }
  }

  if(bjet_pt.second>30)passed=true;

  return passed;


}

bool DoubleHiggsAnalysis::PhotonSelection(TClonesArray *branchPhoton){

  bool passed=false;
  cout<<"branch"<<branchPhoton->GetEntries()<<endl;
  if(branchPhoton->GetEntries() < 2) return passed;
  cout<<"passed branch"<<endl;
  Photon *pho1, *pho2;

  // Take first two photons
  pho1 = (Photon *) branchPhoton->At(0);
  pho2 = (Photon *) branchPhoton->At(1);

  if (pho1->PT < 40. || pho2->PT < 25. )return passed;
  if (pho1->Eta > 2.5 || pho2->Eta > 2.5 )return passed;
  passed=true;

  return passed;

}
