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



double cleverRange(TH1* h,float fac=1.2, float minY=1.e-3)
{
   float maxY =  fac * h->GetBinContent(h->GetMaximumBin());
   //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
   h->SetAxisRange(minY,maxY,"Y");
   return maxY;
}

void handsomeTH1( TH1 *a=0, int col =1, float size=1, int markerstyle=20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();
}


void easyLeg( TLegend *a=0, TString head="") 
{
  a->SetBorderSize(0);
  a->SetHeader(head);
  a->SetTextFont(42);
  //  a->SetTextSize(17);
  a->SetLineColor(1);
  a->SetLineStyle(1);
  a->SetLineWidth(1);
  a->SetFillColor(0);
  a->SetFillStyle(0);

}


void analyzeHistograms() { 
  
  //  TFile* f= new TFile("treeFiles/v5_histograms/trackHistograms_all.root");
  TFile* f= new TFile("treeFiles/v6_histograms/trackHistograms_minNhit3_all.root");

  TH1F* hslopeAXY = (TH1F*)f->Get("hSlopeAXY");
  TH1F* hslopeAZY = (TH1F*)f->Get("hSlopeAZY");
  TH1F* hslopeBYZ = (TH1F*)f->Get("hSlopeBYZ");
  TH1F* hslopeBXZ = (TH1F*)f->Get("hSlopeBXZ");
  TH2D* hax_by = (TH2D*)f->Get("ax_by");
  TH1D* haxResTot = (TH1D*)f->Get("axResTot");

  TH1D* haxRes[10];
  for ( int yid=2 ; yid<6 ; yid++)
    haxRes[yid] = (TH1D*)f->Get(Form("axRes_%d",yid));

  TH2D* htime_z = (TH2D*)f->Get("htime_z");

  TF1 *fun = hslopeAXY->GetFunction("gaus");
  hslopeAXY->GetListOfFunctions()->Remove(fun);
  fun = hslopeAZY->GetFunction("gaus");
  hslopeAZY->GetListOfFunctions()->Remove(fun);
  fun = hslopeBYZ->GetFunction("gaus");
  hslopeBYZ->GetListOfFunctions()->Remove(fun);
  fun = hslopeBXZ->GetFunction("gaus");
  hslopeBXZ->GetListOfFunctions()->Remove(fun);
  fun = haxResTot->GetFunction("gaus");
  haxResTot->GetListOfFunctions()->Remove(fun);


  TH2D* htime_z_run[10];
  TFile* f_run[10];
  for ( int irun = 1 ; irun<=9 ; irun++) {
    if (irun == 4 ) continue;
    f_run[irun] = new TFile(Form("treeFiles/v6_histograms/trackHistograms_minNhit3_run%d.root",irun));
    htime_z_run[irun] = (TH2D*)f_run[irun]->Get("htime_z");
    htime_z_run[irun]->SetName(Form("htime_z_run%d",irun));
  }
  
  
  TCanvas* cvs1 = new TCanvas("cvs1","",800,400);
  cvs1->Divide(2,1);
  cvs1->cd(1);
  
  hslopeAXY->SetAxisRange(-.4,.4,"X");
  cleverRange(hslopeAXY,1.8);
  hslopeAXY->Draw("pe");
  handsomeTH1(hslopeBYZ,2);
  hslopeBYZ->Draw("same pe");
  TLegend *l1 = new TLegend(0.2179929,0.6869624,0.779292,0.9055645,NULL,"brNDC");
  easyLeg(l1,"Slope");
  l1->AddEntry(hslopeBXZ,"BDC","l");
  l1->AddEntry(hslopeAZY,"ATTPC","l");
  l1->Draw();
  
  cvs1->cd(2);
  cleverRange(hslopeAZY,2.5);
  hslopeAZY->SetAxisRange(-.4,.4,"X");
  hslopeAZY->Draw("pe");
  handsomeTH1(hslopeBXZ,2);
  hslopeBXZ->Draw("same pe");
  l1 = new TLegend(0.2179929,0.7869624,0.779292,0.9055645,NULL,"brNDC");
  easyLeg(l1,"Slope");
  l1->AddEntry(hslopeBXZ,"BDC","l");
  l1->AddEntry(hslopeAZY,"ATTPC (v_{drift}=4.6cm/#mus)","l");
  l1->Draw();

  TCanvas* cvs2 = new TCanvas("cvs2","",800,400);
  cvs2->Divide(2,1);
  cvs2->cd(1);
  hax_by->Draw("colz");
  cvs2->cd(2);
  haxResTot->Rebin(4);
  haxResTot->Fit("gaus");
  haxResTot->Draw();


  TH1F* hxRes_vs_yid = new TH1F("hxRes_vs_yid",";ATTPC y index;#Delta (x_{ATTPC} - x_{BDC}) (mm)",8,-.5,7.5);
  TH1F* hxShift_vs_yid = (TH1F*)hxRes_vs_yid->Clone("hxShift_vs_yid");
  TCanvas* cvs3 = new TCanvas("cvs2","",800,400);
  cvs3->Divide(4,2);
  for ( int yid=2 ; yid<6 ; yid++) {
    cvs3->cd(yid+1);
    haxRes[yid]->Rebin(4);
    haxRes[yid]->SetAxisRange(-10,10,"X");
    haxRes[yid]->Draw("pe");
    haxRes[yid]->Fit("gaus");
    TF1 *gausFit = (TF1*) haxRes[yid]->GetFunction("gaus");
    float gausMean = gausFit->GetParameter(1);
    float gausMeanErr = gausFit->GetParError(1);
    float gausSig = gausFit->GetParameter(2);
    float gausSigErr = gausFit->GetParError(2);
    hxShift_vs_yid->SetBinContent(yid+1, gausMean);
    hxShift_vs_yid->SetBinError(yid+1, gausMeanErr);
    hxRes_vs_yid->SetBinContent(yid+1, gausSig);
    hxRes_vs_yid->SetBinError(yid+1, gausSigErr);
  }
  
  TCanvas* cvs4 = new TCanvas("cvs4","",500,500);
  handsomeTH1(hxShift_vs_yid,1);
  handsomeTH1(hxRes_vs_yid,2);
  hxShift_vs_yid->SetAxisRange(-4,4,"Y");
  hxShift_vs_yid->Draw();
  hxRes_vs_yid->Draw("same");
  l1 = new TLegend(0.3393574,0.1957895,1,0.4210526,NULL,"brNDC");
  easyLeg(l1,"x_{ATTPC} - x_{BDC}");
  l1->AddEntry(hxRes_vs_yid,"#sigma","l");
  l1->AddEntry(hxShift_vs_yid,"Mean","l");
  l1->Draw();

  TCanvas* cvs5 = new TCanvas("cvs5","",500,500);
  htime_z->Draw("colz");



  TCanvas* cvs6 = new TCanvas("cvs6","",900,900);
  cvs6->Divide(3,3);
  for ( int irun = 1 ; irun<=9 ; irun++) {
    if (irun == 4 ) continue;
    cvs6->cd(irun);
    htime_z_run[irun]->Draw("colz");
    l1 = new TLegend(0.5943349,0.1346505,1,0.3607903,NULL,"brNDC");
    easyLeg(l1,Form("Run %d",irun));
    l1->Draw();
  }

  TCanvas* cvs7 = new TCanvas("cvs7","",800,400);
  cvs7->Divide(2,1);
  cvs7->cd(1);
  TH2D* htime_z_type1 = (TH2D*)htime_z->Clone("htime_z_type1");
  htime_z_type1->Reset();
  TH2D* htime_z_type2 = (TH2D*)htime_z_type1->Clone("htime_z_type2");
  
  for ( int irun = 1 ; irun<=9 ; irun++) {
    if (irun<=7) htime_z_type1->Add(htime_z_run[irun]);
    else  htime_z_type2->Add(htime_z_run[irun]);
  }
  htime_z_type1->Draw("colz");
  l1 = new TLegend(0.5943349,0.1346505,1,0.3607903,NULL,"brNDC");
  easyLeg(l1,"Run 1 - 7");
  l1->Draw();
    
  cvs7->cd(2);
  htime_z_type2->Draw("colz");
  l1 = new TLegend(0.5943349,0.1346505,1,0.3607903,NULL,"brNDC");
  easyLeg(l1,"Run 8 - 9");
  l1->Draw();
  
  
}


 
