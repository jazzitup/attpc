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

vector<float> doCluster( TH2F* hAdc, TH2F* hTime, float seedThr, TH1F* timediff);

bool isDebugMode = true ;

void treeToTrack( int numEvents = -1 ) {  // # of events to be analyzed.  If -1, we analyze everything


  float seedThr = 100;
  
  TFile* fileIn = new TFile("./trees/treeOfHits_50evts.root");
  TTree* t = (TTree*)fileIn->Get("hit");
  int evtId;
  int nhits;
  float xTree[256];  // x = row * 100mm / 8
  float yTree[256];  // y = col * 100mm / 32
  int xId[256];  // x = row * 100mm / 8
  int yId[256];  // y = col * 100mm / 32
  float timeTree[256];
  float adcTree[256];
  
  t->SetBranchAddress("eventId",&evtId);
  t->SetBranchAddress("nhits",&nhits);
  t->SetBranchAddress("x",xTree);
  t->SetBranchAddress("y",yTree);
  t->SetBranchAddress("xid",xId);
  t->SetBranchAddress("yid",yId);
  t->SetBranchAddress("time",timeTree);
  t->SetBranchAddress("adc",adcTree);
  
  TH1F* hNhits = new TH1F("hNhits",";# of hits;",200,.5,199.5);

  TH2F* hAdcGrid = new TH2F("hAdcGrid_",";xid;yid;ADc",32,-.5,31.5,   8,-.5,7.5);
  TH2F* hTimeGrid = new TH2F("hTimeGrid",";xid;yid;Time",32,-.5,31.5, 8,-.5,7.5);

  TH1F* timediff = new TH1F("timediff","; #Deltat (#mus);",500,-10,10);
  
  auto cvs1 = new TCanvas("cvs1", "", 800, 400);  
  cvs1->Divide(2,1);
  
  int nEvents = t->GetEntries();
  for ( int iev = 0 ; iev <nEvents ; iev++) {
    t->GetEntry(iev);
    hNhits->Fill(nhits);


    hTimeGrid->Reset();
    hAdcGrid->Reset();
    for ( int ihit = 0 ; ihit <nhits ; ihit++)  {
      if ( timeTree[ihit] < 0 ) continue;
      hTimeGrid->SetBinContent( xId[ihit], yId[ihit], timeTree[ihit] * 0.020) ; // micro seconds 
      hAdcGrid->SetBinContent( xId[ihit], yId[ihit], adcTree[ihit] ) ; // micro seconds
      doCluster( hAdcGrid, hTimeGrid, seedThr, timediff);
    }
    if ( isDebugMode) {
      cvs1->cd(1);
      hTimeGrid->Draw("colz");
      cvs1->cd(2);
      hAdcGrid->Draw("colz");
      cvs1->SaveAs(Form("tracking/figure1_%05d.png",evtId));
    }
  }

  auto cvs2 = new TCanvas("cvs2", "", 500, 500);
  if ( isDebugMode ) {
    timediff->Draw();
  }
  
}



vector<float> doCluster( TH2F* hAdc, TH2F* hTime, float seedThr, TH1F* timediff) {
  for ( int iy=0 ; iy<32 ; iy++) {

    // Find the seed in this y strip
    float maxAdc = 0 ; // Max adc in the same-iy-strip
    int maxIx = -1;
    for ( int ix=0 ; ix<32 ; ix++) {
      if  (hAdc->GetBinContent(ix+1, iy+1) > maxAdc)  {
	maxAdc = hAdc->GetBinContent(ix+1, iy+1);
	maxIx = ix;
      }
    }
    if ( maxAdc >  seedThr ) {

      if ( isDebugMode) { 
	if ( maxIx != 0 )
	  timediff->Fill (  hTime->GetBinContent( maxIx, iy+1) - hTime->GetBinContent( maxIx+1, iy+1) );
	if ( maxIx != 7 )
	  timediff->Fill (  hTime->GetBinContent( maxIx+2, iy+1) - hTime->GetBinContent( maxIx+1, iy+1) );
      }
      
      
    }
  }
  vector<float> ret = {0,1} ;
  
  return ret;
  
  
}
