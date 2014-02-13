#include "DoubleHiggsAnalysis.h"

void DoubleHiggsAnalysis::setOutFile(const char *outputFile){

  outFile_=TFile::Open(outputFile,"recreate");

}

void DoubleHiggsAnalysis::setXsec(float xsec){
  xSec_=xsec;
}

void DoubleHiggsAnalysis::setGenEvents(float genevents){
  totalGenEvents_=genevents;
}

void DoubleHiggsAnalysis::setEventWeight(){
  //normalizing to one picobarn
  eventWeight_=xSec_/totalGenEvents_;
}

void DoubleHiggsAnalysis::Analyze(){

  Long64_t numberOfEntries = treeReader_->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader_->UseBranch("Jet");
  TClonesArray *branchPhoton = treeReader_->UseBranch("Photon");
  TClonesArray *branchElectron = treeReader_->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader_->UseBranch("Muon");

  TH1F *histPhoMass = new TH1F("mass", "M_{inv}(#gamma#gamma)", 100, 100.0, 180.0);

  //setXsec
  setXsec(87.);
  setGenEvents(numberOfEntries);
  eventWeight_t=eventWeigh;

  cout<<"starting loop on "<<numberOfEntries<<" entries"<<endl;
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      if(entry%300 == 0)cout<<"### Processing entry: "<<entry<<endl; 
    // Load selected branches with data from specified event
      treeReader_->ReadEntry(entry);
      bool passedJetSelection=false;
      bool passedPhotonSelection=false;
      bool passedLeptonSelection=false;      

      bool looseBtagWP=false;

      Photon *pho1, *pho2;
      // Filling photon variables
      pho1 = (Photon *) branchPhoton->At(0);
      pho2 = (Photon *) branchPhoton->At(1);
      
      passedPhotonSelection=PhotonSelection(branchPhoton);
      if (passedPhotonSelection){
	passedJetSelection=JetSelection(branchJet, looseBtagWP,pho1,pho2);
	passedLeptonSelection=LeptonSelection(branchElectron,branchMuon,pho1,pho2);
      }
      
      //filling histograms and tree
      if(passedJetSelection && passedPhotonSelection){


	ptPhot1_t=pho1->PT;
	ptPhot2_t=pho2->PT;
	etaPhot1_t=pho1->Eta;
	etaPhot2_t=pho2->Eta;
	phiPhot1_t=pho1->Phi;
	phiPhot2_t=pho2->Phi;
	mgg_t=((pho1->P4()) + (pho2->P4())).M();
	histPhoMass->Fill(((pho1->P4()) + (pho2->P4())).M());

	//Filling jet variables
	jet1 = (Jet*) branchJet->At(bjet_indexes.first);
	jet2 = (Jet*) branchJet->At(bjet_indexes.second);

	ptJet_t[0]=jet1->PT;
	ptJet_t[1]=jet2->PT;
	etaJet_t[0]=jet1->Eta;
	etaJet_t[1]=jet2->Eta;
	phiJet_t[0]=jet1->Phi;
	phiJet_t[1]=jet2->Phi;

	deltaRGammaGamma_t=(pho1->P4().DeltaR(pho2->P4()));
	deltaRBB_t=(jet1->P4().DeltaR(jet2->P4()));
	minDeltaRGammaB_t=(pho1->P4().DeltaR(jet1->P4()));
	if(pho1->P4().DeltaR(jet2->P4()) < minDeltaRGammaB_t )minDeltaRGammaB_t=pho1->P4().DeltaR(jet2->P4());
	if(pho2->P4().DeltaR(jet1->P4()) < minDeltaRGammaB_t )minDeltaRGammaB_t=pho2->P4().DeltaR(jet1->P4());
	if(pho2->P4().DeltaR(jet2->P4()) < minDeltaRGammaB_t )minDeltaRGammaB_t=pho2->P4().DeltaR(jet2->P4());
	ptGG_t=(pho1->P4()+pho2->P4()).Pt();
	ptBB_t=(jet1->P4()+jet->P4()).Pt();
	ptBBGG_t=ptBB_t+ptGG_t;
	mgg_t=(pho1->P4()+pho2->P4()).M();
	mgg_t=(jet1->P4()+jet2->P4()).M();
	tree_passedEvents->Fill();
      }
    }
  
  // Show resulting histograms
  //  histJetPT->Write();
  //  histPhoMass->Write();
  //  histPhoLeadPT->Write();
  //  histPhoSubleadPT->Write();
  histPhoMass->Write();
  tree_passedEvents->Write();
  outFile_->Write();
  outFile_->Close();

}

DoubleHiggsAnalysis::DoubleHiggsAnalysis(const char *inputFile,const char *outputFile)
{


  // Create chain of root trees
  chain_=new TChain("Delphes");
  chain_->Add(inputFile);
  
  // Create object of class ExRootTreeReader
  treeReader_ = new ExRootTreeReader(chain_);
  setOutFile(outputFile);  
  createOutputTree();

}


//-------------------------------------------
bool DoubleHiggsAnalysis::JetSelection(TClonesArray *branchJet, bool looseBtagWP,Photon *pho1, Photon *pho2){

  bool passed=false;

  // If event contains at least 2 jet
  if(branchJet->GetEntries() < 2) return passed;
  
  bjet_indexes.first=-1;
  bjet_indexes.second=-1;
  bjet_pt.first=-1.;
  bjet_pt.second=-1.;

  njets_t=0;
  nbjets_t=0;  

  for(int ijet=0;ijet<branchJet->GetEntries();ijet++){
    Jet *jet = (Jet*) branchJet->At(ijet);
    if(jet->P4().DeltaR(pho1->P4())<0.5) continue;
    if(TMath::Abs(jet->Eta)>2.4) continue;
    if(TMath::Abs(jet->PT)<30.) continue;

    bool isBtaggedNormal = (jet->BTag & (1 << 0));
    bool isBtaggedLoose= (jet->BTag & (1 << 1));

    bool isBtagged=isBtaggedNormal;
    if(looseBtagWP)isBtagged=isBtaggedLoose;
    njets_t++;    

    if(!isBtagged) continue;
    nbjets_t++;

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


bool DoubleHiggsAnalysis::LeptonSelection(TClonesArray *branchElectron,TClonesArray *branchMuon,Photon *pho1, Photon *pho2){

  bool isLeptonic=false;
  nleptons_t=0;
  nelectrons_t=0;
  nmuons_t=0;
  minDeltaRElePho_t=-1.;

  if(branchElectron->GetEntries() < 1 || branchMuon->GetEntries() < 1) return isLeptonic;


  electron_index_t=-1;
  ptEle_t=-1;
  for(int ielectron=0;ielectron<branchElectron->GetEntries();ielectron++){
    Electron *electron = (Electron*) branchElectron->At(ielectron);
    if(electron->PT<10.) continue;
    if(electron->Eta<2.4) continue;
    nelectrons_t++;
    if(electron->PT > ptEle_t){
      ptEle_t=electron->PT;
      electron_index_t=ielectron;
    }
  }
  
  if(electron_index_t>-1){
    isLeptonic=1;
    Electron *electron = (Electron*) branchElectron->At(electron_index_t);
    minDeltaRElePho_t=(electron->P4()).DeltaR(pho1->P4());
    if(((electron->P4()).DeltaR(pho2->P4())) < minDeltaRElePho_t) minDeltaRElePho_t=((electron->P4()).DeltaR(pho2->P4()));
    cout<<"#####DeltaR: "<<minDeltaRElePho_t<<endl;
  }

  muon_index_t=-1;
  ptMuon_t=-1;
  for(int imuon=0;imuon<branchMuon->GetEntries();imuon++){
    Muon *muon = (Muon*) branchMuon->At(imuon);
    if(muon->PT<10.) continue;
    if(muon->Eta<2.4) continue;
    nmuons_t++;
    if(muon->PT > ptEle_t){
      ptMuon_t=muon->PT;
      muon_index_t=imuon;
    }
  }
  
  if(muon_index_t>-1){
    isLeptonic=1;
  }

  return isLeptonic;

}

bool DoubleHiggsAnalysis::PhotonSelection(TClonesArray *branchPhoton){

  bool passed=false;
  if(branchPhoton->GetEntries() < 2) return passed;
  Photon *pho1, *pho2;

  // Take first two photons
  pho1 = (Photon *) branchPhoton->At(0);
  pho2 = (Photon *) branchPhoton->At(1);

  if (pho1->PT < 40. || pho2->PT < 25. )return passed;
  if (pho1->Eta > 2.5 || pho2->Eta > 2.5 )return passed;
  passed=true;

  return passed;

}

void DoubleHiggsAnalysis::createOutputTree(){

  tree_passedEvents = new TTree();
  tree_passedEvents->SetName("tree_passedEvents");
  tree_passedEvents->Branch( "eventWeight", &eventWeight_t, "eventWeight_t/F" );
  tree_passedEvents->Branch( "njets", &njets_t, "njets_t/I" );
  tree_passedEvents->Branch( "nbjets", &nbjets_t, "nbjets_t/I" );
  tree_passedEvents->Branch( "ptPhot1", &ptPhot1_t, "ptPhot1_t/F" );
  tree_passedEvents->Branch( "ptPhot2", &ptPhot2_t, "ptPhot2_t/F" );
  tree_passedEvents->Branch( "etaPhot1", &etaPhot1_t, "etaPhot1_t/F" );
  tree_passedEvents->Branch( "etaPhot2", &etaPhot2_t, "etaPhot2_t/F" );
  tree_passedEvents->Branch( "phiPhot1", &phiPhot1_t, "phiPhot1_t/F" );
  tree_passedEvents->Branch( "phiPhot2", &phiPhot2_t, "phiPhot2_t/F" );
  tree_passedEvents->Branch( "mgg", &mgg_t, "mgg_t/F" );
  tree_passedEvents->Branch( "ptJet", ptJet_t, "ptJet_t[njets_t]/F" );
  tree_passedEvents->Branch( "etaJet", etaJet_t, "etaJet_t[njets_t]/F" );
  tree_passedEvents->Branch( "phiJet", phiJet_t, "phiJet_t[njets_t]/F" );
  tree_passedEvents->Branch( "nleptons", &nleptons_t, "nleptons_t/F" );
  tree_passedEvents->Branch( "nelectrons", &nelectrons_t, "nelectrons_t/F" );
  tree_passedEvents->Branch( "nmuons", &nmuons_t, "nmuons_t/F" );
  tree_passedEvents->Branch( "ptEle", &ptEle_t, "ptEle_t/F" );
  tree_passedEvents->Branch( "ptMuon", &ptMuon_t, "ptMuon_t/F" );
  tree_passedEvents->Branch( "ptGG", &ptGG_t, "ptGG_t/F" );
  tree_passedEvents->Branch( "ptBB", &ptBB_t, "ptBB_t/F" );
  tree_passedEvents->Branch( "ptBBGG", &ptBBGG_t, "ptBBGG_t/F" );
  tree_passedEvents->Branch( "minDeltaRElePho", &minDeltaRElePho_t, "minDeltaRElePho_t/F" );
  tree_passedEvents->Branch( "deltaRGammaGamma", &deltaRGammaGamma_t, "deltaRGammaGamma_t/F" );
  tree_passedEvents->Branch( "minDeltaRGammaB", &minDeltaRGammaB_t, "minDeltaRGammaB_t/F" );
  tree_passedEvents->Branch( "deltaRBB", &deltaRBB_t, "deltaRBB_t/F" );
}
