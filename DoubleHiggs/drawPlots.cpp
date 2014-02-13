#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "DrawBase.C"
#include "fitTools.C"

int main(int argc, char* argv[]) {

  if(  argc != 2  ) {
    std::cout << "USAGE: ./drawPlots inputFile" << std::endl;
    exit(23);
  }

  DrawBase* db_nostack = new DrawBase("nostack");

  std::string outputdir_str = "DoubleHiggsPlots/";
  db_nostack->set_outputdir(outputdir_str);
  db_nostack->set_lumiNormalization(3000000.);

  int signalFillColor = 46;
  std::string inputDir="./";

  std::string DoubleHiggsFileName = inputDir + argv[1];
  TFile* DoubleHiggsFile = TFile::Open(DoubleHiggsFileName.c_str());

  db_nostack->add_mcFile( DoubleHiggsFile, "HH", "HH(\\rightarrow bb#gamma#gamma)",signalFillColor+8);
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptPhot1", "eventWeight", 50, 20., 400, "ptPhot1", "p_T Lead Photon (GeV)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "etaPhot1", "eventWeight", 30, -3., 3., "etaPhot1", "#eta Lead Photon (GeV)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptPhot2", "eventWeight", 50, 20., 400, "ptPhot2", "p_T Sublead Photon (GeV)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "etaPhot2", "eventWeight", 30, -3., 3., "etaPhot2", "#eta Lead Photon (GeV)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "njets", "eventWeight", 10, 0.5, 10.5, "njets", "Number of jets","","Entries",true);
  db_nostack->drawHisto_fromTree("tree_passedEvents", "nbjets", "eventWeight", 10, 0.5, 10.5, "nbjets", "Number of b-jets","","Entries",true);
  db_nostack->drawHisto_fromTree("tree_passedEvents", "deltaRGammaGamma", "eventWeight", 15, 0., 5., "deltaRGammaGamma", "#Delta_{R}(#gamma #gamma)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "minDeltaRGammaB", "eventWeight", 20, 0., 5., "mindeltaRGammaB", "#Delta_{R}(#gamma b)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "minDeltaRElePho", "minDeltaRElePho>0.", 20, 0., 4.5, "mindeltaRElePho", "#Delta_{R}(#gamma ele)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "deltaRBB", "eventWeight", 15, 0., 5., "deltaRBB", "#Delta_{R}(bb)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptGG", "eventWeight", 25, 0., 500., "ptGG", "p_{T#gamma#gamma}");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptBB", "eventWeight", 25, 0., 500., "ptBB", "p_{Tbb}");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptBBGG", "eventWeight", 40, 0., 1000., "ptBBGG", "p_{Tbb}+p_{T#gamma#gamma}");
}
