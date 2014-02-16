#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"

#include "TFile.h"
#include "TClonesArray.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "DoubleHiggsAnalysis.C"

using namespace std;

int main(int argc, char *argv[])
{

  if(argc < 2 || argc > 3)
    {
      cout << " Usage: ./DoubleHiggsAnalysis input_file output_file" << endl;
      return 1;
    }
  
  
  DoubleHiggsAnalysis *Analyzer=new DoubleHiggsAnalysis(argv[1],argv[2]);
  Analyzer->Analyze();
  Analyzer->PrintEfficiencies(argv[2]);

}
