#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
//#include "GETAnalyzer.h"
#include "PadMap.h"

void doCluster( TH2F* hAdc, TH2F* hTime, float seedThr, TH1F* timediff, TH1F* hResultX, TH1F* hResultTime, TH1F* hResultNhit, bool useGaus=false);
bool isNearDeadPad( int idx, int idy);

//float bdcTime_to_Sec ( int bdcTime); 

bool isDebugMode = false ;
//bool isDebugMode = true ;
bool isUsingGaus = false;

float dtCut = 0.25;

bool add_BDC_Info = false;


double bY_to_aX ( double bY); 
double bX_to_aZ ( double bX);
double bZ_to_aY ( double bZ);

double aX_to_bY ( double bY); 
double aZ_to_bX ( double bX);
double aY_to_bZ ( double bZ);

PadMap pMap;
//GETAnalyzer analyzer;

void treeToTrack( int numEvents = -1, int runNumber = 1 ) {  // # of events to be analyzed.  If -1, we analyze everything
  
  
  float seedThr = 50;  
  //  float vDrift = 46 ; // in mm/microsecond  <= This must be updated! 
  float vDrift = 55 ; // in mm/microsecond  <= This must be updated! 
  float scTopY = 407.75;
  float scBottomY = -222.25;
  float chi2Cut  = 1;
  
  TString fname = Form("./treeFiles/v7/treeOfHits_run%d_v7_2021Aug30.root",runNumber);
  TFile* fileIn = new TFile(fname);
  //  TFile* fileIn = new TFile("./treeOfHits_muon_run1.root");

  float bucketInMicSec = 0.010;
  TTree* t = (TTree*)fileIn->Get("treeOfHits");

   // Declaration of leaf types
   Int_t           evtId;
   Bool_t          isSpark;
   Double_t        evtTime;
   Double_t        diffTime;
   Int_t           nhits;
   Int_t           run;  //[nhits]
   Double_t        xTree[256];   //[nhits]
   Double_t        yTree[256];   //[nhits]
   Int_t           xId[256];   //[nhits]
   Int_t           yId[256];   //[nhits]
   Int_t           agetId[256];   //[nhits]
   Int_t           chanId[256];   //[nhits]
   Double_t        timeTree[256];   //[nhits]
   Double_t        adcTree[256];   //[nhits]

   // List of branches
   TBranch        *b_evtId;   //!
   TBranch        *b_isSpark;   //!
   TBranch        *b_evtTime;   //!
   TBranch        *b_diffTime;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_run;   //!
   TBranch        *b_xTree;   //!
   TBranch        *b_yTree;   //!
   TBranch        *b_xId;   //!
   TBranch        *b_yId;   //!
   TBranch        *b_timeTree;   //!
   TBranch        *b_adcTree;   //!
   
   /*  int evtId;
  int nhits;
  double evtTime;
  float xTree[256];  // x = row * 100mm / 8
  float yTree[256];  // y = col * 100mm / 32
  int xId[256];  // x = row * 100mm / 8
  int yId[256];  // y = col * 100mm / 32
  float timeTree[256];
  float adcTree[256];
   */
   
  t->SetBranchAddress("eventId", &evtId, &b_evtId);
   t->SetBranchAddress("isSpark", &isSpark, &b_isSpark);
   t->SetBranchAddress("eventTime", &evtTime, &b_evtTime);
   t->SetBranchAddress("diffTime", &diffTime, &b_diffTime);
   t->SetBranchAddress("nHits", &nhits, &b_nhits);
   t->SetBranchAddress("run", &run, &b_run);
   t->SetBranchAddress("x", xTree, &b_xTree);
   t->SetBranchAddress("y", yTree, &b_yTree);
   t->SetBranchAddress("xId", xId, &b_xId);
   t->SetBranchAddress("yId", yId, &b_yId);
   t->SetBranchAddress("time", timeTree, &b_timeTree);
   t->SetBranchAddress("adc", adcTree, &b_adcTree);

   /*   t->SetBranchAddress("eventId",&evtId);
  t->SetBranchAddress("nhits",&nhits);
  t->SetBranchAddress("x",xTree);
  t->SetBranchAddress("y",yTree);
  t->SetBranchAddress("xid",xId);
  t->SetBranchAddress("yid",yId);
  t->SetBranchAddress("time",timeTree);
  t->SetBranchAddress("adc",adcTree);
  t->SetBranchAddress("evtTime", &evtTime); 
*/

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
  Double_t        Zgrad_X[10];   //[trkNumX]
  Double_t        Zgrad_Y[10];   //[trkNumY]
  Double_t        Zc_X[10];   //[trkNumX]
  Double_t        Zc_Y[10];   //[trkNumY]
   
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
   TBranch        *b_Zgrad_X;   //!
   TBranch        *b_Zgrad_Y;   //!
   TBranch        *b_Zc_X;   //!
   TBranch        *b_Zc_Y;   //!
   tBdc->SetBranchAddress("Zgrad_X", Zgrad_X, &b_Zgrad_X);
   tBdc->SetBranchAddress("Zgrad_Y", Zgrad_Y, &b_Zgrad_Y);
   tBdc->SetBranchAddress("Zc_X", Zc_X, &b_Zc_X);
   tBdc->SetBranchAddress("Zc_Y", Zc_Y, &b_Zc_Y);
   TBranch        *b_EvtTime_bdc;   //!
   TBranch        *b_dur_sec;   //!
   TBranch        *b_dur_secDif;   //!
   
   double BDC_evt0_time = 0; 
   if (add_BDC_Info)  {
     if (runNumber == 1)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210810_v7.root");
     if (runNumber == 2)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_530_selecttrack_20210810_v7.root");
     if (runNumber == 3)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_531_selecttrack_20210810_v7.root");
     if (runNumber == 4)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_532_selecttrack_20210810_v7.root");
     if (runNumber == 5)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_533_selecttrack_20210810_v7.root");
     if (runNumber == 6)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_534_selecttrack_20210810_v7.root");
     if (runNumber == 7)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_550_selecttrack_20210810_v7.root");
     if (runNumber == 8)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_560_selecttrack_20210810_v7.root");
     if (runNumber == 9)       tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_561_selecttrack_20210810_v7.root");
     
     
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

   TH1F* hScTopX = new TH1F("hScTopX",";SC top x (mm);",800,-200,200); 
   TH1F* hScBottomX = (TH1F*)hScTopX->Clone("hScBottomX");

   TH1F* hScTopZ = new TH1F("hScTopZ",";SC top z (mm);",800,00,400); 
   TH1F* hScBottomZ = (TH1F*)hScTopZ->Clone("hScBottomZ");

   
   TH2F* hAdc = new TH2F("hAdc",";x (mm);y (mm);ADC",   32, -.5*3.125, 31.5*3.125,   8, -.5*12.5, 7.5*12.5);
   TH2F* hTime = new TH2F("hTime",";x (mm);y (mm);Time",32, -.5*3.125, 31.5*3.125,   8, -.5*12.5, 7.5*12.5);
   
   TH2F* hAdcGrid = new TH2F("hAdcGrid",";xid;yid;ADC",   32,-.5,31.5,   8,-.5,7.5);
   TH2F* hTimeGrid = new TH2F("hTimeGrid",";xid;yid;Time",32,-.5,31.5, 8,-.5,7.5);
   
  TH1F* timediff = new TH1F("timediff","; #Deltat (#mus);",500,-10,10);
  
  TH1F* hResultX = new TH1F("hist_chg"," ; yid ; x_mean",8 , -.5, 7.5);
  TH1F* hResultTime = new TH1F("hist_time"," ; yid ; Timing (#mus)",8 , -.5, 7.5);
  TH1F* hResultNhit = new TH1F("hist_nhits"," ; yid ; N_{hit}",8 , -.5, 7.5);

  TH1F* hSlopeAXY = new TH1F("hSlopeAXY","; dx/dy (ATTPC coor.);",200,-1,1); 
  TH1F* hSlopeAZY = (TH1F*)hSlopeAXY->Clone("hSlopeAZY");
  hSlopeAZY->SetXTitle("dz/dy (ATTPC coor.)");

  TH1F* hSlopeBYZ = (TH1F*)hSlopeAXY->Clone("hSlopeBYZ");
  hSlopeBYZ->SetXTitle("dy/dz (BDC coor.)");
  TH1F* hSlopeBXZ = (TH1F*)hSlopeAXY->Clone("hSlopeBXZ");
  hSlopeBXZ->SetXTitle("dx/dz (BDC coor.)");

  TH1F* hChi2YX = new TH1F("hChi2yx",";#Chi^2 for XY track fit;",100,0,10);
  
  
  TGraph* gResultXY ;
  TGraph* gResultYX ;//  x와 y를 바꾼 것.   피팅할 때는 이것이 더 편함. 
  TGraph* gResultTime ;
  TGraph* gResultTimeYX ;//  x와 y를 바꾼 것.   피팅할 때는 이것이 더 편함. 
  TF1 * fLinYX = new TF1("fLinYX", "[0]+[1]*x", 0, 100); // as a function of Y
  TF1 * fLinTime = new TF1("fLinTime", "[0]+[1]*x", 0, 100); // as a function of Y
  
  
  TH1F* hATTPCTime = new TH1F("hattpctime",";AT-TPC time (s);", 1000,0,100000);
  TH1F* hBDCTime = (TH1F*)hATTPCTime->Clone("hbdctime");
  TH1F* hTimeDiff = new TH1F("hattpctimediff",";t_{AT-TPC} - t_{BDC} (s);", 2000,-1,1);

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
      cout << "minDiff = " << minDiff << endl;
      hTimeDiff->Fill( minDiff);
      //      cout << " attpc_time = " << atime_arr[iev] << "    ";
      //      cout << " bdc_time = " << btime_arr[ index_attpc_to_bdc[iev] ]  <<"    diff = " << minDiff  << endl;
    }
  }
  
  
  
  
  TH2D* ax_by = new TH2D("ax_by",";x_{ATTPC coor.} (mm) ;y_{BDC coor.} (mm)",200,0,100,200,0,100);
  TH2D* adxdy_bdydz = new TH2D("adxdy_bdydz",";dx/dy_{ATTPC coor.} ; -dy/dz_{BDC coor.}",120,-1,1,120,-1,1);
  TH1D* angularRes_adxdy_bdydz = new TH1D("angularREs_adxdy_bdydz",";dx/dy_{ATTPC coor.} - [-dy/dz_{BDC coor.}]",100,-.1,.1);
  TH1D* axResTot = new TH1D("axResTot","; x_{AT-TPC} - x_{BDC} (mm)",200,-20,20);
  TH1D* axRes[10];
  for ( int idy = 0 ; idy<8 ; idy++) {
    axRes[idy] = (TH1D*)axResTot->Clone(Form("axRes_%d",idy));
  }
  TH1D* hxShift_y = new TH1D("hxShift_y","; y (index); <#Delta x> (mm)",8,-.5,7.5);
  TH1D* hxRes_y = new TH1D("hxRes_y","; y (index); #Sigma(#Delta x) (mm)",8,-.5,7.5);
  
  TH2D* htime_z = new TH2D("htime_z",";TPC hit time (#mus); z (mm) from BDC ref.  ",50,0,7,50,20,160);
  TH1D *htimeProf;
  TH1D* hNhitPerCluster = new TH1D("hNhitPerCluster",";N^{hit} per clusters",6,-.5,5.5);
  
  
  
  int nEvents = t->GetEntries();
  
  //  for ( int iev = 450 ; iev <nEvents ; iev++) {
  int attpcCounter=0;
  int bdcCounter=0;
  
  for ( int iev = 0 ; iev <nEvents ; iev++) {
    t->GetEntry(iev);
    hNhits->Fill(nhits);
    if (add_BDC_Info)
      tBdc->GetEntry(index_attpc_to_bdc[iev]);
    
    if (isSpark)
      continue;
    
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
      hTime->    SetBinContent(  xId[ihit]+1, yId[ihit]+1, timeTree[ihit] * bucketInMicSec) ; // micro seconds
      hTimeGrid->SetBinContent(  xId[ihit]+1, yId[ihit]+1,     timeTree[ihit] * bucketInMicSec) ; // micro seconds 
      
    }
    
    doCluster( hAdcGrid, hTimeGrid, seedThr, timediff, hResultX, hResultTime,hResultNhit, isUsingGaus);
    // Fill hNhitPerCluster
    for ( int idy=0; idy<8 ; idy++) { 
      hNhitPerCluster->Fill (hResultNhit->GetBinContent(idy+1) );
    }
    
    
    float px[8];
    float py[8];
    float ptime[8];
    int nClus=0;
    for ( int iy=2 ; iy<6 ; iy++) {
      //      if ( hResultX->GetBinContent(iy+1) > 0 ) {
      if ( hResultNhit->GetBinContent(iy+1) >= 3 ) { // at least three hits
	px[nClus] = hResultX->GetBinContent(iy+1) * 3.125; 
	py[nClus] = iy * 12.5 ;
	ptime[nClus] = hResultTime->GetBinContent(iy+1);
	
	if (add_BDC_Info) {
	  if (index_attpc_to_bdc[iev] == -1)
	    cout << " BDC not triggered" << endl;
	  else {
	    //	  if (trckNumX==1) {
	    if ( (trckNumY==1) && (trckNumX==1)) {
	      // The aY position is  py[nClus]
	      
	      double bZ = aY_to_bZ(py[nClus]);
	      
	      //      bdc_z = Ygrad[0] * bdc_y + Yc[0] ;
	      //        // z = Ygrad[ix] *y + Yc ;
	      double bY = (bZ - Yc[0])/Ygrad[0];
	      double bX = (bZ - Xc[0])/Xgrad[0];
	      double aX = bY_to_aX(bY);
	      double aZ =  bX_to_aZ(bX);
	      
	      //	    if ( px[nClus] >= 7*3.125 && px[nClus] < 16*3.125 )  {
	      //    if ( hResultNhit->GetBinContent(idy+1) >=3 ) {
	      if ( (px[nClus] >= 20 && px[nClus] < 90) )   {
		ax_by->Fill( aX, px[nClus]);
		axResTot->Fill ( px[nClus] - aX );
		axRes[iy]->Fill ( px[nClus] - aX );
		htime_z->Fill ( ptime[nClus], aZ);
	      }
	    }
	    
	  }
	} //if (add_BDC_Info)	
	
	nClus++; 	
      }
    }
    
    
    if ( nClus <3 )
      continue;   // We don't need to fit this!!
    
    attpcCounter++;
    if (trckNumY==1) 
      bdcCounter++;
    
    gResultXY =  new TGraph(nClus,px,py);
    gResultYX =  new TGraph(nClus,py,px);
    
    fLinYX->SetParameters(  0 ,    0);
    //    fLinYX->SetParLimits(0,1,100);
    fLinYX->SetParLimits(1,-1,1);
    gResultYX->Fit(fLinYX, "M R Q");

    double yxChi2 = fLinYX->GetChisquare() ; 
    hChi2YX->Fill(yxChi2);


    gResultTime =  new TGraph(nClus,ptime,py);
    gResultTimeYX =  new TGraph(nClus,py,ptime);
    
    fLinTime->SetParameters(4,0);
    gResultTimeYX->Fit(fLinTime, "M R Q");

    float theX0   = fLinYX->GetParameter(0);
    float theDxdy = fLinYX->GetParameter(1);
    float scTopX = theX0 + scTopY*theDxdy;
    float scBottomX = theX0 + scBottomY*theDxdy;

    float theT0 = fLinTime->GetParameter(0);
    float theDtdy = fLinTime->GetParameter(1);
    float scTopT = theT0 + scTopY*theDtdy;
    float scBottomT = theT0 + scBottomY*theDtdy;
        
    if (yxChi2 < chi2Cut)  {
      hScTopX->Fill (scTopX);
      hScBottomX->Fill (scBottomX);
      hScTopZ->Fill (scTopT * vDrift);  // ns * mm/ns
      hScBottomZ->Fill (scBottomT * vDrift);
    }
    // angle :
    //if ( (trckNumY==1) && (trckNumX==1)) {
    hSlopeAXY->Fill(fLinYX->GetParameter(1));
    hSlopeAZY->Fill ( fLinTime->GetParameter(1)*vDrift );
    
    if (add_BDC_Info) {
      if (trckNumY==1) {
	hSlopeBYZ->Fill( -Zgrad_X[0]);
	hSlopeBXZ->Fill( Zgrad_Y[0]);
	adxdy_bdydz->Fill( fLinYX->GetParameter(1),  -Zgrad_X[0]);
	angularRes_adxdy_bdydz->Fill ( fLinYX->GetParameter(1) + Zgrad_X[0]);
      }
    }


    
    if ( isDebugMode) {
      TCanvas* cvs1 = new TCanvas("cvs1", "", 800, 800);  
      cvs1->Divide(2,2);
      cvs1->cd(2);
      hTime->Draw("colz");
      cvs1->cd(1);
      hAdc->Draw("colz");
      gResultXY->SetMarkerSize(1);
      gResultXY->SetMarkerColor(kWhite);
      gResultXY->Draw("same p");
      //                aY = - 410.25 +  Ygrad[0] * (aX + 48.5625)  + Yc[0]
      TF1 *f1 = new TF1("f1","-410.25 + [0]*(x + 48.5625)  + [1]",0,100);
      if ( add_BDC_Info && (trckNumY==1)) {
	f1->SetParameter(0,  Ygrad[0]);
	f1->SetParameter(1,  Yc[0]);
	f1->SetLineColor(4);
	f1->SetLineWidth(4);
	f1->Draw("same");
      }
      
      cvs1->cd(3);
      //      hResultX->SetAxisRange(0,32,"Y");
      //      hResultX->Draw("e");
      htemp->Reset();
      htemp->SetXTitle("y (mm)");
      htemp->SetYTitle("x (mm)");
      htemp->SetAxisRange(0,100,"Y");
      htemp->DrawCopy();
      gResultYX->Draw("same p");
      cvs1->cd(4);
      htemp->Reset();
      htemp->SetYTitle("Time");
      htemp->SetAxisRange(0,10,"Y");
      htemp->DrawCopy();
      gResultTimeYX->Draw("same p");
      
      cvs1->SaveAs(Form("tracking/v7/figure1_%05d_run%d_BDCevent%d.png",evtId,runNumber,(int)Event));
      delete f1;
    }
  }
  

  //  if ( isDebugMode ) {
  TCanvas* cvs0 = new TCanvas("cvs0", "", 500,500);
  hChi2YX->Draw();
  //  }
  
  if ( isDebugMode && add_BDC_Info) {
    TCanvas* cvs2 = new TCanvas("cvs22", "", 500, 500);
    timediff->Draw();
  }
  
  
  TCanvas* cvs3 = new TCanvas("cvs3", "", 800, 800);
  cvs3->Divide(2,2);
  cvs3->cd(1);
  hSlopeAXY->Fit("gaus");
  hSlopeAXY->SetStats(1);
  hSlopeAXY->Draw();
  cvs3->cd(2);
  hSlopeAZY->Fit("gaus");
  hSlopeAZY->SetStats(1);
  hSlopeAZY->Draw();
  TLegend *l1 = new TLegend(0.2179929,0.7869624,0.779292,0.9455645,NULL,"brNDC");
  l1->SetHeader(Form("Assumed v_{drift} = %.1f cm/#mus",vDrift*0.1));
  l1->SetBorderSize(0);
  l1->SetLineColor(1);
  l1->SetLineStyle(1);
  l1->Draw();
  
  cvs3->cd(3);
  hSlopeBYZ->Fit("gaus");
  hSlopeBYZ->SetStats(1);
  hSlopeBYZ->Draw();
  cvs3->cd(4);
  hSlopeBXZ->Fit("gaus");
  hSlopeBXZ->SetStats(1);
  hSlopeBXZ->Draw();
  
  
  /*  TCanvas* cvs4 = new TCanvas("cvs4","",800,400);
  cvs4->Divide(2,1);
  cvs4->cd(1);
  ax_by->Draw("colz");
  cvs4->cd(2);
  axResTot->Draw();
  axResTot->Fit("gaus");
  TF1 *theGausFit = (TF1*) axResTot->GetFunction("gaus");
  float gausMean = theGausFit->GetParameter(1);
  float gausSig = theGausFit->GetParameter(2);
  TLegend *l3 = new TLegend(0.2179929,0.2869624,0.879292,0.5455645,NULL,"brNDC");
  l3->SetHeader("Fit result");
  l3->AddEntry("",Form("Mean = %.1f mm", gausMean),"");
  l3->AddEntry("",Form("#sigma = %.1f mm", gausSig),"");
  l3->SetBorderSize(0);
  l3->SetLineColor(1);
  l3->SetLineStyle(1);
  l3->Draw(); 
  */
  TCanvas* cvs45 = new TCanvas("cvs45","",1200,600);
  if (add_BDC_Info) { cvs45->Divide(4,2);
    for ( int iy = 0 ; iy<8 ; iy++) {
      cvs45->cd(iy+1);
      axRes[iy]->Draw();
      //    axRes[iy]->Fit("gaus");
      //    TF1 *gausFit = (TF1*) axRes[iy]->GetFunction("gaus");
      //  float gausMean = gausFit->GetParameter(1);
    //    float gausMeanErr = gausFit->GetParError(1);
    //    float gausSig = gausFit->GetParameter(2);
    //    float gausSigErr = gausFit->GetParError(2);
    //    hxShift_y->SetBinContent(iy+1, gausMean);
    //    hxShift_y->SetBinError(iy+1, gausMeanErr);
    //    hxRes_y->SetBinContent(iy+1, gausSig);
    //    hxRes_y->SetBinError(iy+1, gausSigErr);
    //    
    /*    TLegend *l31 = new TLegend(0.2179929,0.2869624,0.879292,0.5455645,NULL,"brNDC");
    l31->SetHeader("Fit result");
    l31->AddEntry("",Form("Mean = %.1f mm", gausMean),"");
    l31->AddEntry("",Form("#sigma = %.1f mm", gausSig),"");
    l31->SetBorderSize(0);
    l31->SetLineColor(1);
    l31->SetLineStyle(1);
    l31->Draw(); */
  }
  }
  
  TCanvas* cvs46 = new TCanvas("cvs46","",500,500);
  TLegend *l31 = new TLegend(0.2179929,0.2869624,0.879292,0.5455645,NULL,"brNDC");
  
  if (add_BDC_Info) {
    hxShift_y->SetMarkerColor(1);
    hxShift_y->SetLineColor(1);
    hxRes_y->SetMarkerColor(2);
    hxRes_y->SetLineColor(2);
    
    hxShift_y->Draw();
    hxRes_y->Draw("same");
    l31->SetHeader("t_{ATTPC} - t_{BDC}");
    l31->AddEntry("hxShift_y","mean","pe");
    l31->AddEntry("hxRes_y","resolution","pe");
    l31->SetBorderSize(0);
    l31->SetLineColor(1);
    //  l31->SetLineStyle(1);
    l31->Draw();
  }
  
  if (add_BDC_Info) {
    TCanvas* cvs5 = new TCanvas("cvs5","",800,400);
    cvs5->Divide(2,1);
    cvs5->cd(1);
    htime_z->Draw("colz");
    cvs5->cd(2);
    htimeProf = (TProfile*)htime_z->ProfileX()->ProjectionX();
    htimeProf->Draw();
    htimeProf->Fit("pol1","M","",1,3);
    
    TF1 *thePol1Fit = (TF1*) htimeProf->GetFunction("pol1");
    float zeroPoint = thePol1Fit->GetParameter(0);
    float drftVel = thePol1Fit->GetParameter(1);
    
    TLegend *l2 = new TLegend(0.2179929,0.2869624,0.879292,0.5455645,NULL,"brNDC");
    l2->SetHeader("Fit result");
    l2->AddEntry("",Form("v_{drift} = %.1f cm/#mus", drftVel*0.1),"");
    l2->AddEntry("",Form("Time_{z=0} = %.1f #mus", zeroPoint*0.1),"");
    l2->SetBorderSize(0);
    l2->SetLineColor(1);
    l2->SetLineStyle(1);
    l2->Draw();

  TCanvas* cvsTime = new TCanvas("cvstime","",800,400);
  cvsTime->Divide(2,1);
  cvsTime->cd(1);
  hATTPCTime->SetLineColor(2);
  hATTPCTime->SetMarkerColor(2);
  hBDCTime->Draw();
  hATTPCTime->Draw("same p");
  cvsTime->cd(2);
  hTimeDiff->Draw();
  }

  
  TFile* fout = new TFile(Form("treeFiles/v7_histograms/trackHistograms_minNhit3_run%d.root",runNumber),"recreate");
  hNhitPerCluster->Write();
  hSlopeAXY->Write();
  hSlopeAZY->Write();
  hScTopX->Write();
  hScBottomX->Write();
  hScTopZ->Write();
  hScBottomZ->Write();
  
  if (add_BDC_Info) {
    hSlopeBYZ->Write();
    hSlopeBXZ->Write();
    ax_by->Write();
    adxdy_bdydz->Write();
    angularRes_adxdy_bdydz->Write();
    axResTot->Write();
    {
      for ( int iy = 0 ; iy<8 ; iy++) {
	axRes[iy]->Write();
      }
    }
    hxShift_y->Write();
    hxRes_y->Write();
    htime_z->Write();
    htimeProf->Write();
  }
  
  fout->Close();

  cout << "* Run No. " << runNumber << endl; 
  cout << "    Good AT-TPC events     : " << attpcCounter << endl;
  cout << "    (Out of it) BDC events : " << bdcCounter << endl;
}



void doCluster( TH2F* hAdc, TH2F* hTime, float seedThr, TH1F* timediff, TH1F* hResultX, TH1F* hResultTime, TH1F* hResultNhit, bool useGaus) {

  hResultX->Reset();
  hResultTime->Reset();
  hResultNhit->Reset();
  TH1F* hForFit = new TH1F("hForFit","", 32,-0.5,31.5);
  TF1* fGaus = new TF1("fGaus","gaus",0,32);
  
  for ( int iy=0 ; iy<8 ; iy++) {
    
    // Find the seed in this y strip
    float maxAdc = 0 ; // Max adc in the same-iy-strip
    float maxAdcTime = -100;
    int maxIx = -1;
    for ( int ix=0 ; ix<32 ; ix++) {
      if  (hAdc->GetBinContent(ix+1, iy+1) > maxAdc)  {
	maxAdc = hAdc->GetBinContent(ix+1, iy+1);
	maxIx = ix;
      }
    }
    if ( maxAdc < seedThr )
      continue;
    if  ( isNearDeadPad( maxIx, iy) )   // if it's near the dead channel
      continue; 
    
    maxAdcTime =  hTime->GetBinContent( maxIx+1, iy+1) ; // Max adc in the same-iy-strip
    if ( isDebugMode) { 
      if ( maxIx != 0 )
	timediff->Fill (  hTime->GetBinContent( maxIx, iy+1) - maxAdcTime );
      if ( maxIx != 7 )
	timediff->Fill (  hTime->GetBinContent( maxIx+2, iy+1) - maxAdcTime ) ;
    }
    
    
    int nValHit = 0;
    float sumAdc = 0 ;
    hForFit->Reset();
    for ( int iclu = -2 ; iclu <=2 ; iclu++) {
      int cluIx = maxIx + iclu;
      if ( (cluIx>-1) && (cluIx<32) && (fabs(hTime->GetBinContent(cluIx+1, iy+1) - maxAdcTime) < dtCut) ) { // cluIx : 0 ~ 31
	if (hAdc->GetBinContent ( cluIx+1, iy+1) > seedThr)
	  {
	    nValHit++;
	    sumAdc = sumAdc + hAdc->GetBinContent ( cluIx+1, iy+1) ;
	    hForFit->SetBinContent( cluIx+1,  hAdc->GetBinContent ( cluIx+1, iy+1) );  
	  }
      }
    }
    //      cout << "fGaus->GetParameter(0)" << fGaus->GetParameter(0) << endl; 
    float meanX = hForFit->GetMean();
    fGaus->SetParameter(0, maxAdc);
    fGaus->SetParameter(1, meanX);
    
    hForFit->Fit(fGaus, "LL M Q Q R");
    
    float finalX =  meanX;
    if (useGaus) finalX = fGaus->GetParameter(1) ;

    //if ( nValHit == 1 )
    //      finalX = maxIx;
    
    hResultX->SetBinContent( iy+1, finalX);
    hResultX->SetBinError( iy+1, 0.001);
    hResultTime->SetBinContent(iy+1, maxAdcTime);
    hResultTime->SetBinError(iy+1, 0.001);
    hResultNhit->SetBinContent(iy+1, nValHit);
  }
  
  delete fGaus;
  delete hForFit;
  //  delete h2dtemp;
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
bool isNearDeadPad( int idx, int idy) {
  bool result = false; 
  //  cout << "idx, idy = " << idx << ", " << idy<<endl;
  if (pMap.IsDeadPad( idx, idy))
    result = true;

  if (idx != 0) {
    if (pMap.IsDeadPad(idx-1,idy)) {
      result = true;
    }
  }

  if (idx != 31) {
    if (pMap.IsDeadPad(idx+1,idy))
      {
	result = true;
      }
  }
  return result;
}

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
//   bdc_z = Ygrad[0] * bdc_y + Yc[0]
//  aY = - 410.25 +  Ygrad[0] * (aX + 48.5625)  + Yc[0] 
