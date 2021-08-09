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

void doCluster( TH2F* hAdc, TH2F* hTime, float seedThr, TH1F* timediff, TH1F* hResChg, TH1F* hResTime);

//float bdcTime_to_Sec ( int bdcTime); 

//bool isDebugMode = true ;
bool isDebugMode = false ;

float dtCut = 0.25;

bool add_BDC_Info = true;


void treeToTrack( int numEvents = -1 ) {  // # of events to be analyzed.  If -1, we analyze everything
  
  
  float seedThr = 100;
  
  float  vDrift = 48 ; // in mm/microsecond  <= This must be updated! 
  TFile* fileIn = new TFile("./treeOfHits_muon_run1.root");

  TTree* t = (TTree*)fileIn->Get("hit");
  int evtId;
  int nhits;
  double evtTime;
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
  t->SetBranchAddress("evtTime", &evtTime);

  t->GetEntry(0);  // 0점 맞추기
  double ATTPC_evt0_time = evtTime; 

  
  TChain* tBdc = new TChain("trkTree");
  Int_t           Event;
  Int_t           trckNumX;
  Int_t           trckNumY;
  Double_t        Xgrad[10];   //[trckNumX]
  Double_t        Ygrad[10];   //[trckNumY]
  Double_t        Xc[10];   //[trckNumX]
  Double_t        Yc[10];   //[trckNumY]
  //  Int_t           EvtTime_bdc;
  //  Double_t        dur_sec;
  Double_t        dur_secDif;
   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_trckNumX;   //!
   TBranch        *b_trckNumY;   //!
   TBranch        *b_Xgrad;   //!
   TBranch        *b_Ygrad;   //!
   TBranch        *b_Xc;   //!
   TBranch        *b_Yc;   //!
   TBranch        *b_EvtTime_bdc;   //!
   TBranch        *b_dur_sec;   //!
   TBranch        *b_dur_secDif;   //!
   
   double BDC_evt0_time = 0; 
   if (add_BDC_Info)  {
     tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210808_v3.root");
     tBdc->SetBranchAddress("Event", &Event, &b_Event);
     tBdc->SetBranchAddress("trckNumX", &trckNumX, &b_trckNumX);
     tBdc->SetBranchAddress("trckNumY", &trckNumY, &b_trckNumY);
     tBdc->SetBranchAddress("Xgrad", Xgrad, &b_Xgrad);
     tBdc->SetBranchAddress("Ygrad", Ygrad, &b_Ygrad);
     tBdc->SetBranchAddress("Xc", Xc, &b_Xc);
     tBdc->SetBranchAddress("Yc", Yc, &b_Yc);
     //     tBdc->SetBranchAddress("EvtTime", &EvtTime_bdc, &b_EvtTime_bdc);
     //     tBdc->SetBranchAddress("dur_sec", &dur_sec, &b_dur_sec); 
     tBdc->SetBranchAddress("dur_secDif", &dur_secDif, &b_dur_secDif);
	
     tBdc->GetEntry(0);
     cout << " BDC reference time (str) = " << dur_secDif << endl;
     cout << " BDC reference time (sec) = " << BDC_evt0_time << endl; 
   }
   
   
   
   
   
   TH2F* htemp = new TH2F("htemp",";Y (mm); x (mm)",100,0,100, 100,0,100);
   TH1F* hNhits = new TH1F("hNhits",";# of hits;",200,.5,199.5);
  
  TH2F* hAdc = new TH2F("hAdc",";x (mm);y (mm);ADC",   32, -.5*3.125, 31.5*3.125,   8, -.5*12.5, 7.5*12.5);
  TH2F* hTime = new TH2F("hTime",";x (mm);y (mm);Time",32, -.5*3.125, 31.5*3.125,   8, -.5*12.5, 7.5*12.5);

  TH2F* hAdcGrid = new TH2F("hAdcGrid",";xid;yid;ADC",   32,-.5,31.5,   8,-.5,7.5);
  TH2F* hTimeGrid = new TH2F("hTimeGrid",";xid;yid;Time",32,-.5,31.5, 8,-.5,7.5);

  TH1F* timediff = new TH1F("timediff","; #Deltat (#mus);",500,-10,10);

  TH1F* hResultCharge = new TH1F("hist_chg"," ; yid ; x_mean",8 , -.5, 7.5);
  TH1F* hResultTime = new TH1F("hist_time"," ; yid ; Timing (#mus)",8 , -.5, 7.5);

  TH1F* hThetaXY = new TH1F("hThetaXY","; Angle in XY plane (degree);",50,-50,50); 
  TH1F* hThetaYZ = (TH1F*)hThetaXY->Clone("hThetaYZ");
  
  TGraph* gResultXY ;
  TGraph* gResultYX ;//  x와 y를 바꾼 것.   피팅할 때는 이것이 더 편함. 
  TGraph* gResultTime ;
  TGraph* gResultTimeYX ;//  x와 y를 바꾼 것.   피팅할 때는 이것이 더 편함. 
  TF1 * fLinYX = new TF1("fLinYX", "[0]+[1]*x", 0, 100); // as a function of Y
  TF1 * fLinTime = new TF1("fLinTime", "[0]+[1]*x", 0, 100); // as a function of Y
  
  
  TH1F* hATTPCTime = new TH1F("hattpctime",";AT-TPC time (s);", 1000,0,100000);
  TH1F* hBDCTime = (TH1F*)hATTPCTime->Clone("hbdctime");

  const int maxEvents = 10000;
  int index_attpc_to_bdc[maxEvents];
  for ( int i=0; i<maxEvents; i++)   index_attpc_to_bdc[i] = -1;
  
  if (add_BDC_Info)  { 
    cout << "BDC timing" << endl;
    for ( int jev =0  ; jev < tBdc->GetEntries() ; jev++) {
      tBdc->GetEntry(jev);
      double bdc_time  = dur_secDif - BDC_evt0_time ;
      //       cout << " EvtTime_bdc = " << bdc_time << endl;
      hBDCTime->Fill ( bdc_time);
    }
  }
  
  for ( int iev = 0 ; iev <t->GetEntries() ; iev++) {
    t->GetEntry(iev);
    double attpc_time = evtTime - ATTPC_evt0_time ;
    hATTPCTime->Fill ( attpc_time);     
  }
  
  

   
   
   TCanvas* cvsTime = new TCanvas("cvstime","",400,400);
   hATTPCTime->SetLineColor(2);
   hATTPCTime->SetMarkerColor(2);
   hBDCTime->Draw();
   hATTPCTime->Draw("same p");
   
     return;
   
   
   int nEvents = t->GetEntries();
   //  for ( int iev = 450 ; iev <452 ; iev++) {
   //  for ( int iev = 450 ; iev <nEvents ; iev++) {
   for ( int iev = 0 ; iev <nEvents ; iev++) {
     t->GetEntry(iev);
    hNhits->Fill(nhits);

    double attpc_time = evtTime - ATTPC_evt0_time ; 
    cout << " AT-TPC event time = " << attpc_time << endl; 
    
     

    hTimeGrid->Reset();
    hAdcGrid->Reset();
    hTime->Reset();
    hAdc->Reset();
    for ( int ihit = 0 ; ihit <nhits ; ihit++)  {
      if ( timeTree[ihit] < 0 ) continue;


      hAdc->    SetBinContent(  xId[ihit]+1, yId[ihit]+1, adcTree[ihit] ) ; // micro seconds
      hAdcGrid->SetBinContent( xId[ihit]+1, yId[ihit]+1,     adcTree[ihit] ) ; // micro seconds

      hTime->    SetBinContent(  xId[ihit]+1, yId[ihit]+1, timeTree[ihit] * 0.020) ; // micro seconds
      hTimeGrid->SetBinContent(  xId[ihit]+1, yId[ihit]+1,     timeTree[ihit] * 0.020) ; // micro seconds 

      
      
    }
    
    doCluster( hAdcGrid, hTimeGrid, seedThr, timediff, hResultCharge, hResultTime);
    
    int nClus=0;
    float px[8];
    float py[8];
    float ptime[8];
    for ( int iy=0 ; iy<8 ; iy++) {
      if ( hResultCharge->GetBinContent(iy+1) > 0 ) {
	px[nClus] = hResultCharge->GetBinContent(iy+1) * 3.125; 
	py[nClus] = iy * 12.5 ;
	ptime[nClus] = hResultTime->GetBinContent(iy+1);
	nClus++; 
      }
    }
    if ( nClus <3 )
      continue;   // We don't need to fit this!!
  
    gResultXY =  new TGraph(nClus,px,py);
    gResultYX =  new TGraph(nClus,py,px);
    
    fLinYX->SetParameters(  0 ,    0);
    //    fLinYX->SetParLimits(0,1,100);
    fLinYX->SetParLimits(1,-1,1);
    gResultYX->Fit(fLinYX, "M R Q");

    gResultTime =  new TGraph(nClus,ptime,py);
    gResultTimeYX =  new TGraph(nClus,py,ptime);

    fLinTime->SetParameters(4,0);
    gResultTimeYX->Fit(fLinTime, "M R Q");

    // angle :
    float thetaXY = atan(fLinYX->GetParameter(1)) * 180 / 3.141592; 
    hThetaXY->Fill(thetaXY);

    float thetaYZ = atan(fLinTime->GetParameter(1)*vDrift) * 180 / 3.141592;  
    hThetaYZ->Fill(thetaYZ);

    
    if ( isDebugMode) {
      TCanvas* cvs1 = new TCanvas("cvs1", "", 800, 800);  
      cvs1->Divide(2,2);
      cvs1->cd(1);
      hTime->Draw("colz");
      cvs1->cd(2);
      hAdc->Draw("colz");
      gResultXY->SetMarkerSize(1);
      gResultXY->SetMarkerColor(kWhite);
      gResultXY->Draw("same p");
      cvs1->cd(3);
      //      hResultCharge->SetAxisRange(0,32,"Y");
      //      hResultCharge->Draw("e");
      htemp->Reset();
      htemp->SetXTitle("y(mm)");
      htemp->SetYTitle("Charge");
      htemp->SetAxisRange(0,100,"Y");
      htemp->DrawCopy();
      gResultYX->Draw("same p");
      cvs1->cd(4);
      htemp->Reset();
      htemp->SetYTitle("Time");
      htemp->SetAxisRange(0,10,"Y");
      htemp->DrawCopy();
      gResultTimeYX->Draw("same p");

      cvs1->SaveAs(Form("tracking/figure1_%05d.png",evtId));
    }

  }
  
  if ( isDebugMode ) {
    TCanvas* cvs2 = new TCanvas("cvs22", "", 500, 500);
    timediff->Draw();
  }

  TCanvas* cvs3 = new TCanvas("cvs3", "", 800, 400);
  cvs3->Divide(2,1);
  cvs3->cd(1);
  hThetaXY->Fit("gaus");
  hThetaXY->SetStats(1);
  hThetaXY->Draw();
  cvs3->cd(2);
  hThetaYZ->Fit("gaus");
  hThetaYZ->SetStats(1);
  hThetaYZ->Draw();
  cout << "hThetaXY->Integral()=" << hThetaXY->Integral() << endl;
  
}



void doCluster( TH2F* hAdc, TH2F* hTime, float seedThr, TH1F* timediff, TH1F* hResChg, TH1F* hResTime) {

  hResChg->Reset();
  hResTime->Reset();
  TH1F* hForFit = new TH1F("hForFit","", 32,-0.5,31.5);
  TF1* fGaus = new TF1("fGaus","gaus",0,32);

  
  for ( int iy=0 ; iy<8 ; iy++) {
    

    // Find the seed in this y strip
    float maxAdc = 0 ; // Max adc in the same-iy-strip
    int maxIx = -1;
    for ( int ix=0 ; ix<32 ; ix++) {
      if  (hAdc->GetBinContent(ix+1, iy+1) > maxAdc)  {
	maxAdc = hAdc->GetBinContent(ix+1, iy+1);
	maxIx = ix;
      }
    }
    float maxAdcTime =  hTime->GetBinContent( maxIx+1, iy+1) ; // Max adc in the same-iy-strip
  
    if ( maxAdc >  seedThr ) {
      
      if ( isDebugMode) { 
	if ( maxIx != 0 )
	  timediff->Fill (  hTime->GetBinContent( maxIx, iy+1) - maxAdcTime );
	if ( maxIx != 7 )
	  timediff->Fill (  hTime->GetBinContent( maxIx+2, iy+1) - maxAdcTime ) ;
      }
      
      float sumAdc = 0 ;
      float xWgtSumAdc = 0 ;
      float meanX ;

      hForFit->Reset();
      int nValHit = 0;
      for ( int iclu = -2 ; iclu <=2 ; iclu++) {
	int cluIx = maxIx + iclu;
	if ( (cluIx>-1) && (cluIx<32) && (fabs(hTime->GetBinContent(cluIx+1, iy+1) - maxAdcTime) < dtCut) ) { // cluIx : 0 ~ 31
	  sumAdc = sumAdc + hAdc->GetBinContent ( cluIx+1, iy+1) ;
	  xWgtSumAdc = xWgtSumAdc +  hAdc->GetBinContent ( cluIx+1, iy+1) * cluIx ; 
	  
	  hForFit->SetBinContent( cluIx+1,  hAdc->GetBinContent ( cluIx+1, iy+1) );  
	  nValHit++;
	}
      }

      //      cout << "fGaus->GetParameter(0)" << fGaus->GetParameter(0) << endl; 
      meanX = xWgtSumAdc/sumAdc;
      fGaus->SetParameter(0, maxAdc);
      fGaus->SetParameter(0, meanX);
      
      hForFit->Fit(fGaus, "LL M Q Q R");

      //      cout << " nValHit = " << nValHit << endl;
      //       cout << "mean from gaus. fit  = " << fGaus->GetParameter(1) << endl; 
      //      cout << "arithmatic mean      = " << meanX << endl;

      //      cout << "Avg - seed = " << meanX - maxIx << endl;

      float finalX =  fGaus->GetParameter(1) ;
      if ( nValHit == 1 )
	finalX = maxIx;
      //      finalX = meanX;
      
      hResChg->SetBinContent( iy+1, finalX);
      hResChg->SetBinError( iy+1, 0.001);
      hResTime->SetBinContent(iy+1, maxAdcTime);
      hResTime->SetBinError(iy+1, 0.001);
      if ( 1==0)  {
	TCanvas* tempcvs = new TCanvas("tempcvs","",500,500);
	hForFit->Draw();
	tempcvs->SaveAs(Form("tracking/figure1_iy%d.png",iy));
	delete tempcvs;
      }
    }
  }
  
  delete fGaus;
  delete hForFit;
  
}

/* float bdcTime_to_Sec ( int bdcTime) {
  // DD HH MM SS
  int dd = bdcTime/1000000;
  int hh = (bdcTime%1000000) /10000;
  int mm = (bdcTime%10000) / 100;
  int ss = bdcTime%100;

  return dd*86400 + hh*3600 + mm*60 + ss;
  
}
*/
