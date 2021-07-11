#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DataFrame.h"
#include "PadMap.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TH1S.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

void treeToTrack( int numEvents = -1 ) {  // # of events to be analyzed.  If -1, we analyze everything
  
  TFile* fileIn = new TFile("./trees/treeOfHits_1kEvents.root");
  TTree* t = (TTree*)fileIn->Get("hit");
  int evtId;
  int nhits;
  float xTree[256];  // x = row * 100mm / 8
  float yTree[256];  // y = col * 100mm / 32
  float timeTree[256];
  float adcTree[256];
  
  t->SetBranchAddress("eventId",&evtId);
  t->SetBranchAddress("nhits",&nhits);
  t->SetBranchAddress("x",xTree);
  t->SetBranchAddress("y",yTree);
  t->SetBranchAddress("time",timeTree);
  t->SetBranchAddress("adc",adcTree);
  
  TH1F* hTime = new TH1F("hTime",";timing (#mus);",512,0,512*0.020);
  TH1F* hNhits = new TH1F("hNhits",";# of hits;",200,.5,199.5);
  int nEvents = t->GetEntries();
  for ( int iev = 0 ; iev <nEvents ; iev++) {
    t->GetEntry(iev);
    hNhits->Fill(nhits);
    for ( int ihit = 0 ; ihit <nhits ; ihit++) 
      hTime->Fill( timeTree[ihit] * 0.020) ; // micro seconds
  }
    auto cvs1 = new TCanvas("cvs1", "", 800, 400);
    cvs1->Divide(2,1);
    cvs1->cd(1);
    hNhits->Draw();
    cvs1->cd(2);
    hTime->Draw();
}

