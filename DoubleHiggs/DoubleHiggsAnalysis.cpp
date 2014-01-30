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

int main()
{

  cout<<"starting"<<endl;
  DoubleHiggsAnalysis *Analyzer=new DoubleHiggsAnalysis("delphes_out_SM_14TeV_NoPileUp.root");
  cout<<"analyzing"<<endl;
  Analyzer->Analyze();


}
