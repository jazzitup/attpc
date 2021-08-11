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

double bY_to_aX ( double bY);
double bX_to_aZ ( double bX);
double bZ_to_aY ( double bZ);

double aX_to_bY ( double bY);
double aZ_to_bX ( double bX);
double aY_to_bZ ( double bZ);


void bdcAna( int numEvents = -1 ) {  // # of events to be analyzed.  If -1, we analyze everything
  
  
  float seedThr = 100;

  float  vDrift = 48 ; // in mm/microsecond  <= This must be updated! 
  TChain* tBdc = new TChain("trkTree");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_530_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_531_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_532_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_533_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_534_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_550_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_560_selecttrack_20210810_v5.root");
  tBdc->Add("BDCTrackingData/bdcAnaTrack_Data_SJ_Run_561_selecttrack_20210810_v5.root");

   Int_t           Event;
   Int_t           trckNumX;
   Int_t           trckNumY;
   Double_t        Xgrad[1];   //[trkNumX]
   Double_t        Ygrad[1];   //[trkNumY]
   Double_t        Xc[1];   //[trkNumX]
   Double_t        Yc[1];   //[trkNumY]
   Double_t        Zgrad_X[1];   //[trkNumX]
   Double_t        Zgrad_Y[1];   //[trkNumY]
   Double_t        Zc_X[1];   //[trkNumX]
   Double_t        Zc_Y[1];   //[trkNumY]
   Int_t           EvtTime;
   Double_t        dur_sec;
   Int_t           EvtTimeDif;
   Double_t        dur_secDif;
   Int_t           EvtTag;
      
   // List of branches
    TBranch        *b_Event;   //!
   TBranch        *b_trkNumX;   //!
   TBranch        *b_trkNumY;   //!
   TBranch        *b_Xgrad;   //!
   TBranch        *b_Ygrad;   //!
   TBranch        *b_Xc;   //!
   TBranch        *b_Yc;   //!
   TBranch        *b_Zgrad_X;   //!
   TBranch        *b_Zgrad_Y;   //!
   TBranch        *b_Zc_X;   //!
   TBranch        *b_Zc_Y;   //!
   TBranch        *b_EvtTime;   //!
   TBranch        *b_dur_sec;   //!
   TBranch        *b_EvtTimeDif;   //!
   TBranch        *b_dur_secDif;   //!
   TBranch        *b_EvtTag;   //!


   tBdc->SetBranchAddress("Event", &Event, &b_Event);
   tBdc->SetBranchAddress("trckNumX", &trckNumX, &b_trkNumX);
   tBdc->SetBranchAddress("trckNumY", &trckNumY, &b_trkNumY);
   tBdc->SetBranchAddress("Xgrad", Xgrad, &b_Xgrad);
   tBdc->SetBranchAddress("Ygrad", Ygrad, &b_Ygrad);
   tBdc->SetBranchAddress("Xc", Xc, &b_Xc);
   tBdc->SetBranchAddress("Yc", Yc, &b_Yc);
   tBdc->SetBranchAddress("Zgrad_X", Zgrad_X, &b_Zgrad_X);
   tBdc->SetBranchAddress("Zgrad_Y", Zgrad_Y, &b_Zgrad_Y);
   tBdc->SetBranchAddress("Zc_X", Zc_X, &b_Zc_X);
   tBdc->SetBranchAddress("Zc_Y", Zc_Y, &b_Zc_Y);
   tBdc->SetBranchAddress("EvtTime", &EvtTime, &b_EvtTime);
   tBdc->SetBranchAddress("dur_sec", &dur_sec, &b_dur_sec);
   tBdc->SetBranchAddress("EvtTimeDif", &EvtTimeDif, &b_EvtTimeDif);
   tBdc->SetBranchAddress("dur_secDif", &dur_secDif, &b_dur_secDif);
   tBdc->SetBranchAddress("EvtTag", &EvtTag, &b_EvtTag);

      
   //   TH1F* hNhits = new TH1F("hNhits",";# of hits;",200,.5,199.5);

   int nBdcEvents = tBdc->GetEntries();
   cout << " Entries = " << nBdcEvents << endl;
   int counts = 0 ; 
   int countsX = 0 ; 
   int countsY = 0 ;
   //   float zCenter = 614;
   float zValue =  454 ; //454; //aY_to_bZ ( 50) ;//454; //454;
   TH1D* hBx = new TH1D("hBx",Form(";BDC x when z = %d (mm)",(int)zValue),50,0,200);
   TH1D* hBy = new TH1D("hBy",Form(";BDC y when z = %d (mm)",(int)zValue),50,0,200);
   TH1D* hAx = new TH1D("hAx",Form(";ATTPC x when y = %d (mm)",(int)(bZ_to_aY(zValue))),50,0,200);
   TH1D* hAz = new TH1D("hAz",Form(";ATTPC z when y = %d (mm)",(int)(bZ_to_aY(zValue))),50,0,200);
   TH1D* hBXslope = new TH1D("hBXslope",";BDC dx/dz;",50,-.3,.3);
   TH1D* hBYslope = new TH1D("hBYslope",";BDC dy/dz;",50,-.3,.3);
   TH1D* hAXslope = new TH1D("hAXslope",";ATTPC dx/dy;",50,-.3,.3);
   TH1D* hAZslope = new TH1D("hAZslope",";ATTPC dz/dy;",50,-.3,.3);

   TH2D* hx2 = new TH2D("hx2","",20,0,200,20,-.3,.3);
   TH2D* hy2 = new TH2D("hy2","",20,0,200,20,-.3,.3);
      

   cout << "===== Tracks ======"<< endl;
   for ( int jev = 0 ; jev <nBdcEvents ; jev++) {
     tBdc->GetEntry(jev);			  
     if ( (trckNumX == 1)&&(trckNumY==1)) 
       counts++;
     if  (trckNumX == 1) countsX++;
     if  (trckNumY == 1) countsY++;
     
     if ( ! ((trckNumX == 1)&&(trckNumY==1)))
       continue;
       
     // when aY = -100 mm 
     float theAY = -100;
     float theAZ;
     float theAX;
     float velAY =1.0;
     float velAZ;
     float velAX;
     

     // XZ plane
     for ( int ix = 0 ; ix<trckNumX;  ix++) {
       // x = Zgrad_X  * z + Zc_X ; 
       float xValue = Zgrad_X[ix] * zValue + Zc_X[ix];
       //       cout << "xCenter = " << xCenter << endl;
       hBx->Fill( xValue);
       hAz->Fill ( bX_to_aZ (xValue) );
       hBXslope->Fill( Zgrad_X[ix]);
       hAZslope->Fill( -Zgrad_X[ix]);
       hx2->Fill ( xValue, Zgrad_X[ix]);
       
       theAZ = bX_to_aZ(xValue) + (-Zgrad_X[ix]) * ( theAY - bZ_to_aY (zValue) ) ;
       velAZ =  -Zgrad_X[ix] ; 
       
     }
     // YZ plane
     for ( int iy = 0 ; iy<trckNumY;  iy++) {
       // x = Zgrad_X  * z + Zc_X ; 
       float yValue = Zgrad_Y[iy] * zValue + Zc_Y[iy];
       hBy->Fill( yValue);
       hAx->Fill ( bY_to_aX (yValue) );
       hBYslope->Fill( Zgrad_Y[iy]);
       hAXslope->Fill(  Zgrad_Y[iy]);
       hy2->Fill ( yValue, Zgrad_Y[iy]);

       theAX = bY_to_aX(yValue) + (+Zgrad_Y[iy]) * ( theAY - bZ_to_aY (zValue) ) ;
       velAX =  Zgrad_Y[iy];
     }
     
     float velNorm = sqrt (velAX*velAX + velAY*velAY + velAZ*velAZ);
     velAX = velAX/velNorm;
     velAY = velAY/velNorm;
     velAZ = velAZ/velNorm;
     
     cout << theAX <<" " << theAY << " " << theAZ << "  " ;
     cout << velAX <<" " << velAY << " " << velAZ << "  " << endl;
     
     // print out x,y,z when z = -1
     //    double bY_to_aX ( double bY) {
     //double bX_to_aZ ( double bX) {
     //double bZ_to_aY ( double bZ) {
     //     float theAZ = -10;
     //     float ax = bY_to_aX ( Zgrad_Y[iy] * theZ + Zc_Y[iy] ) ;
     //    float ay = bZ_to_aY ( Zgrad_Y[iy] * theZ + Zc_Y[iy] ) ;
     
     //     cout << " 
     
   }
			
			
   TCanvas* cvs1 = new TCanvas("cvs1", "", 800, 800);
   cvs1->Divide(2,2);
   cvs1->cd(1);
   hBx->Draw();
   hBx->Fit("gaus");
   cvs1->cd(2);
   hBy->Draw();
   hBy->Fit("gaus");
   cvs1->cd(3);

   hBXslope->Draw();
   hBXslope->Fit("gaus");
   cvs1->cd(4);

   hBYslope->Draw();
   hBYslope->Fit("gaus");
   cvs1->Update();
			
   TCanvas* cvs2 = new TCanvas("cvs2", "", 800, 800);
   cvs2->Divide(2,2);
   cvs2->cd(1);
   hAz->Draw();
   hAz->Fit("gaus");
   cvs2->cd(2);
   hAx->Draw();
   hAx->Fit("gaus");
   cvs2->cd(3);

   hAZslope->Draw();
   hAZslope->Fit("gaus");
   cvs2->cd(4);
   hAXslope->Draw();
   hAXslope->Fit("gaus");
   cvs2->Update();

			
			
   TCanvas* cvs3 = new TCanvas("cvs3", "", 800, 400);
   cvs3->Divide(2,1);
   cvs3->cd(1);
   hx2->Draw("colz");
   cvs3->cd(2);
   hy2->Draw("colz");
   
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
