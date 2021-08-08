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

bool isDebugMode = true ;
//bool isDebugMode = false ;

void bdcAna( int numEvents = -1 ) {  // # of events to be analyzed.  If -1, we analyze everything


  float seedThr = 100;

  float  vDrift = 48 ; // in mm/microsecond  <= This must be updated! 
  TChain* tBdc = new TChain("trkTree");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_520_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_530_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_531_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_532_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_533_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_534_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_550_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_560_20210806_v2.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_561_20210806_v2.root");
  Int_t           Event;
  Int_t           trckNumX;
  Int_t           trckNumY;
  Double_t        Xgrad[10];   //[trckNumX]
  Double_t        Ygrad[10];   //[trckNumY]
  Double_t        Xc[10];   //[trckNumX]
  Double_t        Yc[10];   //[trckNumY]
   Int_t           EvtTime;
   Double_t        dur_sec;
   
   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_trckNumX;   //!
   TBranch        *b_trckNumY;   //!
   TBranch        *b_Xgrad;   //!
   TBranch        *b_Ygrad;   //!
   TBranch        *b_Xc;   //!
   TBranch        *b_Yc;   //!
   TBranch        *b_EvtTime;   //!
   TBranch        *b_dur_sec;   //!
   
   tBdc->SetBranchAddress("Event", &Event, &b_Event);
   tBdc->SetBranchAddress("trckNumX", &trckNumX, &b_trckNumX);
   tBdc->SetBranchAddress("trckNumY", &trckNumY, &b_trckNumY);
   tBdc->SetBranchAddress("Xgrad", Xgrad, &b_Xgrad);
   tBdc->SetBranchAddress("Ygrad", Ygrad, &b_Ygrad);
   tBdc->SetBranchAddress("Xc", Xc, &b_Xc);
   tBdc->SetBranchAddress("Yc", Yc, &b_Yc);
   tBdc->SetBranchAddress("EvtTime", &EvtTime, &b_EvtTime);
   tBdc->SetBranchAddress("dur_sec", &dur_sec, &b_dur_sec);
   
   //   TH1F* hNhits = new TH1F("hNhits",";# of hits;",200,.5,199.5);

   int nBdcEvents = tBdc->GetEntries();
   cout << " Entries = " << nBdcEvents << endl;
   int counts = 0 ; 
   int countsX = 0 ; 
   int countsY = 0 ;
   //   float zCenter = 614;
   float zCenter = 454;
   TH1D* hXCenter = new TH1D("hXCenter",Form(";x when z = %d (mm)",(int)zCenter),40,0,200);
   TH1D* hYCenter = new TH1D("hYCenter",Form(";y when z = %d (mm)",(int)zCenter),40,0,200);
   TH1D* hXgrad = new TH1D("hXgrad",";dx/dz;",50,-.3,.3);
   TH1D* hYgrad = new TH1D("hYgrad",";dy/dz;",50,-.3,.3);
   for ( int jev = 0 ; jev <nBdcEvents ; jev++) {
     tBdc->GetEntry(jev);			  
     if ( (trckNumX == 1)&&(trckNumY==1)) 
       counts++;
     if  (trckNumX == 1) countsX++;
     if  (trckNumY == 1) countsY++;
     

     // xz plane : 
     for ( int ix = 0 ; ix<trckNumX;  ix++) {
       // z = Xgrad[ix] *x + Xc ; 
       float xCenter = (zCenter - Xc[ix]) / Xgrad[ix] ;
       cout << "xCenter = " << xCenter << endl;
       hXCenter->Fill( xCenter);
       hXgrad->Fill( 1/Xgrad[ix]);
     }
     //yz plane:
     for ( int iy = 0 ; iy<trckNumY;  iy++) {
       // z = Ygrad[ix] *y + Yc ; 
       float yCenter = (zCenter - Yc[iy]) / Ygrad[iy] ;
       cout << "yCenter = " << yCenter << endl;
       hYCenter->Fill(yCenter);
       hYgrad->Fill( 1/Ygrad[iy]);
     }
     
     
   }
   cout << " total events = " << nBdcEvents << endl;
   cout << " X fired events = " << countsX << endl;
   cout << " Y fired events = " << countsY << endl;
   cout << " Both X Y fired events = " << counts << endl;

     TCanvas* cvs1 = new TCanvas("cvs1", "", 800, 800);
     cvs1->Divide(2,2);
     cvs1->cd(1);
     hXCenter->Draw();
     cvs1->cd(2);
     hYCenter->Draw();
     cvs1->cd(3);
     hXgrad->Draw();
     cvs1->cd(4);
     hYgrad->Draw();

}

   /*
     
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
   */
