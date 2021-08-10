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


double bY_to_aX ( double bY); 
double bX_to_aZ ( double bX);
double bZ_to_aY ( double bZ);

double aX_to_bY ( double bY); 
double aZ_to_bX ( double bX);
double aY_to_bZ ( double bZ);





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
     tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210810_v5.root");
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
  TH1F* hTimeDiff = new TH1F("hattpctime",";AT-TPC time (s);", 2000,-20,20);

  const int maxEvents = 10000;
  double atime_arr[maxEvents]; // AT-TPC time array
  double btime_arr[maxEvents]; // BDC time array
  int index_attpc_to_bdc[maxEvents];        
  
  if (add_BDC_Info)  { 
    for ( int iev = 0 ; iev <t->GetEntries() ; iev++) {
      t->GetEntry(iev);
      atime_arr[iev] =  evtTime - ATTPC_evt0_time;
    }
    for ( int jev =0  ; jev < tBdc->GetEntries() ; jev++) {
      tBdc->GetEntry(jev);
      btime_arr[jev]  = dur_secDif - BDC_evt0_time ;
    }
  }

  for ( int iev = 0 ; iev <t->GetEntries() ; iev++) {
    double minDiff = 10;
      int matchedIndex =  -1;
      for ( int jev =0  ; jev < tBdc->GetEntries() ; jev++) {
	if ( fabs(atime_arr[iev] - btime_arr[jev]) < fabs(minDiff) ) {
	  minDiff = atime_arr[iev] - btime_arr[jev] ;
	  matchedIndex = jev; 
	}
      }
      index_attpc_to_bdc[iev] = matchedIndex ;
      hTimeDiff->Fill( minDiff);
      //      cout << " attpc_time = " << atime_arr[iev] << "    ";
      //      cout << " bdc_time = " << btime_arr[ index_attpc_to_bdc[iev] ]  <<"    diff = " << minDiff  << endl;
  }
  
  
  
  
  TCanvas* cvsTime = new TCanvas("cvstime","",800,400);
  cvsTime->Divide(2,1);
  cvsTime->cd(1);
  hATTPCTime->SetLineColor(2);
  hATTPCTime->SetMarkerColor(2);
  hBDCTime->Draw();
  hATTPCTime->Draw("same p");
  cvsTime->cd(2);
  hTimeDiff->Draw();
  //  return;

  TH2D* ax_bx = new TH2D("ax_bx",";ax;bx",100,-100,100,100,-100,100);
  TH1D* axRes = new TH1D("axRes","",30,-4,4);
  
  int nEvents = t->GetEntries();
  //  for ( int iev = 450 ; iev <452 ; iev++) {
  //  for ( int iev = 450 ; iev <nEvents ; iev++) {
  for ( int iev = 0 ; iev <nEvents ; iev++) {
    t->GetEntry(iev);
    hNhits->Fill(nhits);
    tBdc->GetEntry(index_attpc_to_bdc[iev]);
    
    double attpc_time = evtTime - ATTPC_evt0_time ; 
    //    cout << " AT-TPC event time = " << attpc_time << endl; 

    
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

	if (index_attpc_to_bdc[iev] == -1)
	  cout << " BDC not triggered" << endl;
	else {
	  if ( (trckNumY==1) && (trckNumX==1)) {
	    // The aY position is  py[nClus]

	    double bZ = aY_to_bZ(py[nClus]);

	    //      bdc_z = Ygrad[0] * bdc_y + Yc[0] ;
	    //        // z = Ygrad[ix] *y + Yc ;
	    double bY = (bZ - Yc[0])/Ygrad[0];
	    double aX = bY_to_aX(bY);
	    //	    cout << " cluster x = " << px[nClus] << ",    BDC ref = " << aX << ",   diff = " << px[nClus]-aX << endl;
	    ax_bx->Fill( aX, px[nClus]);
	    if ( px[nClus] > 40 && px[nClus] < 45 ) 
	      axRes->Fill ( px[nClus] - aX );
	    // double bY_to_aX ( double bY);
	    // double bX_to_aZ ( double bX);
	    // double bZ_to_aY ( double bZ);
	    
	    //      bdc_z = Ygrad[0] * bdc_y + Yc[0] ;
	    //      float bdc_y = (bdc_z - Yc[0]) / Ygrad[0] ;
	  }
	  
	}	
	
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

  TCanvas* cvs4 = new TCanvas("cvs4","",800,400);
  cvs4->Divide(2,1);
  cvs4->cd(1);
  ax_bx->Draw("colz");
  cvs4->cd(2);
  axRes->Draw();
  axRes->Fit("gaus");
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

double bY_to_aX ( double bY) {
  return bY - (97 - 48.4375);   // bdc y = 97 -> attpc x = 48.4375
}                                                                       // before : 270-155.5 = 114.5 | after : 270-173 = 97
double bX_to_aZ ( double bX) {
  return  -bX + 173.5 ;                         // before : 171, after : 171+2.5 = 173.5
}
double bZ_to_aY ( double bZ) {
  return  bZ  - (454 - 43.75) ;
}

double aX_to_bY ( double aX) {
  return aX + (97 - 48.4375);   // bdc y = 97 -> attpc x = 48.4375
}
double aZ_to_bX ( double aZ) {
  return -aZ + 173.5 ;
}
double aY_to_bZ ( double aY) {
  return  aY  + (454 - 43.75) ;
}
