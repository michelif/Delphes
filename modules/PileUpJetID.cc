
/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables
 *
 *  \author S. Zenz
 *
 */

#include "modules/PileUpJetID.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PileUpJetID::PileUpJetID() :
  fItJetInputArray(0),fTrackInputArray(0),fNeutralInputArray(0)
{

}

//------------------------------------------------------------------------------

PileUpJetID::~PileUpJetID()
{

}

//------------------------------------------------------------------------------

void PileUpJetID::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);
  fParameterR = GetDouble("ParameterR", 0.5);
  fUseConstituents = GetInt("UseConstituents", 0);

  // import input array(s)

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();


  //  cout << "BeforE SCZ additions in init" << endl;
  //  cout << GetString("TrackInputArray", "ParticlePropagator/tracks") << endl;
  //  cout << GetString("EFlowTrackInputArray", "ParticlePropagator/tracks") << endl;

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "ParticlePropagator/tracks"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();


  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

  //  cout << " end of INIT " << endl;

}

//------------------------------------------------------------------------------

void PileUpJetID::Finish()
{
  //  cout << "In finish" << endl;

  if(fItJetInputArray) delete fItJetInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;

}

//------------------------------------------------------------------------------

void PileUpJetID::Process()
{
  //  cout << "start of process" << endl;

  Candidate *candidate, *constituent;
  TLorentzVector momentum, area;

  //  cout << "BeforE SCZ additions in process" << endl;

  // SCZ
  Candidate *trk;

  // loop over all input candidates
  fItJetInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    momentum = candidate->Momentum;
    area = candidate->Area;

    float sumpt = 0.;
    float sumptch = 0.;
    float sumptchpv = 0.;
    float sumptchpu = 0.;
    float sumdrsqptsq = 0.;
    float sumptsq = 0.;
    int nc = 0;
    int nn = 0;
    float pt_ann[5];

    for (int i = 0 ; i < 5 ; i++) {
      pt_ann[i] = 0.;
    }

    if (fUseConstituents) {
      TIter itConstituents(candidate->GetCandidates());
      while((constituent = static_cast<Candidate*>(itConstituents.Next()))) {
        float pt = constituent->Momentum.Pt();
        float dr = candidate->Momentum.DeltaR(constituent->Momentum);
	sumpt += pt;
	sumdrsqptsq += dr*dr*pt*pt;
	sumptsq += pt*pt;
	if (constituent->Charge == 0) {
	  nn++;
	} else {
	  if (constituent->IsRecoPU) {
	    sumptchpu += pt;
	  } else {
	    sumptchpv += pt;
	  }
	  sumptch += pt;
	  nc++;
	}
	for (int i = 0 ; i < 5 ; i++) {
	  if (dr > 0.1*i && dr < 0.1*(i+1)) {
	    pt_ann[i] += pt;
	  }
	}
      }
    } else {
      // Not using constituents, using dr
      fItTrackInputArray->Reset();
      while ((trk = static_cast<Candidate*>(fItTrackInputArray->Next()))) {
	if (trk->Momentum.DeltaR(candidate->Momentum) < fParameterR) {
	  float pt = trk->Momentum.Pt();
	  sumpt += pt;
	  sumptch += pt;
	  if (trk->IsRecoPU) {
	    sumptchpu += pt;
	  } else {
	    sumptchpv += pt;
	  }
	  float dr = candidate->Momentum.DeltaR(trk->Momentum);
	  sumdrsqptsq += dr*dr*pt*pt;
	  sumptsq += pt*pt;
	  nc++;
	  for (int i = 0 ; i < 5 ; i++) {
	    if (dr > 0.1*i && dr < 0.1*(i+1)) {
	      pt_ann[i] += pt;
	    }
	  }
	}
      }
      fItNeutralInputArray->Reset();
      while ((trk = static_cast<Candidate*>(fItNeutralInputArray->Next()))) {
	if (trk->Momentum.DeltaR(candidate->Momentum) < fParameterR) {
	  float pt = trk->Momentum.Pt();
	  sumpt += pt;
	  float dr = candidate->Momentum.DeltaR(trk->Momentum);
	  sumdrsqptsq += dr*dr*pt*pt;
	  sumptsq += pt*pt;
	  nn++;
	  for (int i = 0 ; i < 5 ; i++) {
	    if (dr > 0.1*i && dr < 0.1*(i+1)) {
            pt_ann[i] += pt;
	    }
	  }
	}
      }
    }

    if (sumptch > 0.) {
      candidate->Beta = sumptchpu/sumptch;
      candidate->BetaStar = sumptchpv/sumptch;
    } else {
      candidate->Beta = -999.;
      candidate->BetaStar = -999.;
    }
    if (sumptsq > 0.) {
      candidate->MeanSqDeltaR = sumdrsqptsq/sumptsq;
    } else {
      candidate->MeanSqDeltaR = -999.;
    }
    candidate->NCharged = nc;
    candidate->NNeutrals = nn;
    if (sumpt > 0.) {
      candidate->PTD = TMath::Sqrt(sumptsq) / sumpt;
      for (int i = 0 ; i < 5 ; i++) {
        candidate->FracPt[i] = pt_ann[i]/sumpt;
      }
    } else {
      candidate->PTD = -999.;
      for (int i = 0 ; i < 5 ; i++) {
        candidate->FracPt[i] = -999.;
      }
    }

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
