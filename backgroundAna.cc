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
#include "TH1F.h"
#include "TH1S.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

// void DynamicExec();
// double BackgroundFcn(double *x, double *par);
double SignalFcn(double *x, double *par);
// double TotalFcn(double *x, double *par);


using std::cout;
using std::endl;

int OddOrEven( int aget, int chan ) {
  int ret =100;

  //  if ( aget !=2 ) return 100;
  
  if (chan <11)  {
    if (chan%2 ==0)  ret = 0;
    else ret = 1;
  }
  else if ( chan <22 ){
    if (chan%2 ==0)  ret = 1;
    else ret = 0;
  }
  
  else if ( chan <45 ){
    if (chan%2 ==0)  ret = 0;
    else ret = 1;
  }
  else if ( chan <56 ){
    if (chan%2 ==0)  ret = 1;
    else ret = 0;    
  }
  else {
    if (chan%2 ==0)  ret = 0;
    else ret = 1;   
  }      
  
  if ( aget !=0 ) {
    ret = fabs(ret-1);
  }
  return ret;
}

void trimFirstFive(TH1* h ) {
  h->SetBinContent(1,0);
  h->SetBinContent(2,0);
  h->SetBinContent(3,0);
  h->SetBinContent(4,0);
  h->SetBinContent(5,0);
}  
  

void backgroundAna(int startEvent=-1, int endEvent = -1, std::string strInFileName = "RawFiles/cosmic_202104/CoBo_2021-05-13T19h00m11.474s_0000.graw") {

  //  gStyle->SetPalette(kGreyScale);
    
  //void backgroundAna( std::string strInFileName = "RawFiles/cosmic_202104/CoBo_2021-05-13T19h00m11.474s_0000.graw") {
  std::ifstream inFile;
    inFile.open(strInFileName.c_str(), std::ios::in | std::ios::binary);
    if (!inFile.is_open()) {
        cout << "File open error : " << strInFileName << endl;
        exit(-1);
    }
    int threshold1 = 600 ;   // If adc's from all channel are smaller than this, it's considered as the potnetial background.
    int threshold2 = 200 ;   // threshold applied to the hits in the final event display
    auto hPoly = new TH2Poly("Read-out Pad", "", -3.5625, 100.4375, -8.25, 95.75);    // Event display with simple threshold 
    auto hPolyFinal = new TH2Poly("Read-out Pad after background subtraction", "", -3.5625, 100.4375, -8.25, 95.75);    // Event display after proper background subtraction.
    auto c1 = new TCanvas("c1", "c1", 700, 700);
    auto c2 = new TCanvas("c2", "c2", 700, 700);
    auto *hAdcTime = new TH1F("hAdcTime", "", 512, 0, 512);
    TH1F *hADCcollection[32][8];
    TH1F *hADCcollectionBkgSubt[32][8];
    


    int label[32][8]; // signal or bakcgrdoun;
    for ( int i=0;i<32;i++)
      for ( int j=0;j<8;j++)
	label[i][j] = -1;    // -1: undefined, 0: background,   1: signal+background
    
    // auto fitFcn = new TF1("fitFcn", TotalFcn, 0, 512, 4);
    auto fitFcn = new TF1("fitFcn", SignalFcn, 0, 512, 3);

    TH1F *hSig = new TH1F("hSig","",100,0,2000);
    TH1F *hBkg= new TH1F("hBkg","",100,0,2000);

    int nRelAdcBins= 300;
    TH2F *BkgProfileOdd = new TH2F("bkgprofileOdd","",512,0,512,nRelAdcBins,0,nRelAdcBins);
    TH2F *BkgProfileEven = new TH2F("bkgprofileEven","",512,0,512,nRelAdcBins,0,nRelAdcBins);

    TH2F *BkgProfileOddCorr = new TH2F("bkgprofileOddCorr","",512,0,512,nRelAdcBins,0,nRelAdcBins);
    TH2F *BkgProfileEvenCorr = new TH2F("bkgprofileEvenCorr","",512,0,512,nRelAdcBins,0,nRelAdcBins);

    
    DataFrame frame;
    PadMap pMap(1);
    pMap.BuildPad(hPoly);
    pMap.BuildPad(hPolyFinal);

    // Output TTree
    TTree* treeOut = new TTree("hits","a Tree with hits");
    int nhitsForTree;
    float x[256];  // x = row * 100mm / 8
    float y[256];  // y = col * 100mm / 32 
    float time[256];
    float adc[256];
    
    treeOut->Branch("nhits", &nhitsForTree, "nhits/I");
    treeOut->Branch("x", x, "x[nhits]/F");
    treeOut->Branch("y", y, "y[nhits]/F");
    treeOut->Branch("time", time, "time[nhits]/F");
    treeOut->Branch("adc", adc, "adc[nhits]/F");
    
    
    int nEvents = 0;
    while (!inFile.eof()) {
      frame.ReadHeader(inFile);
      if (inFile.eof()) {
	cout << "[EOF]" << endl;
	break;
      }
      frame.ReadItem(inFile);
      
      if ( nEvents > endEvent )
	break;
      
      if ( (startEvent == -1 ) || (nEvents >= startEvent))	  {
	hPoly->ClearBinContents();
	
	int nHits = 0;
	  cout << nEvents << endl;
	  hAdcTime->Reset();
	  BkgProfileEvenCorr->Reset();
	  BkgProfileOddCorr->Reset();
	  BkgProfileOdd->Reset();
	  BkgProfileEven->Reset();
	  
	  for (int colN = 0; colN < 32; colN++) {
	    for (int rowN = 0; rowN < 8; rowN++) {
	      int maxValue, maxBin;
	      for (int buckIdx = 0; buckIdx < 512; buckIdx++) {
		int agetIdx = pMap.GetAgetIdx(colN, rowN);
		int chanIdx = pMap.GetChanIdx(colN, rowN);
		int ADC;
		ADC = frame.GetADC(agetIdx, chanIdx, buckIdx);
		hAdcTime->SetBinContent(buckIdx + 1, ADC);
	      }
	      maxBin = hAdcTime->GetMaximumBin();
	      maxValue = hAdcTime->GetBinContent(maxBin);
	      hADCcollection[colN][rowN] = (TH1F*)hAdcTime->Clone(Form("adc_col%d_row%d",colN,rowN));
	      c1->cd();
	      hAdcTime->SetStats(0);
	      hAdcTime->SetAxisRange(0, 2000, "Y");
	      c1->Update();
	      hAdcTime->Draw();
	      //	      c1->SaveAs(Form("./figures/Event%d_%d_%d.png", nEvents, pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)));
	      if (maxValue > threshold1 ) {
		hPoly->Fill(pMap.GetX(colN), pMap.GetY(colN, rowN), maxValue);
		nHits++;
		
		for (int buckIdx = 0; buckIdx < 512; buckIdx++) {
		  hSig->Fill ( hAdcTime->GetBinContent(buckIdx+1));
		}
		label[colN][rowN] = 1;
		
	      } else {
		hPoly->Fill(pMap.GetX(colN), pMap.GetY(colN, rowN), 1);
		
		for (int buckIdx = 0; buckIdx < 512; buckIdx++) {
		  hBkg->Fill ( hAdcTime->GetBinContent(buckIdx+1));
		}
		label[colN][rowN] = 0;
		
		
		
	      }
            }
	  }
	  c2->cd();
	  c2->SetLogz();
	  hPoly->SetStats(0);
	  hPoly->SetAxisRange(1, 5000, "Z");
	  hPoly->SetTitle(Form("Event #%d;Y [mm];Z [mm]", nEvents));
	  hPoly->Draw("colz");
	  c2->SetRightMargin(0.12);
	  c2->Update();
	  if (nHits > 4) {
	    c2->SaveAs(Form("./figures/track_%d.png", nEvents));
	  }
	
	  auto c3 = new TCanvas("c3", "c3", 500, 500);
	  hSig->SetLineColor(2);
	  hSig->Scale( 1./hSig->Integral());
	  hBkg->Scale( 1./hBkg->Integral());
	  hSig->Draw("hist");
	  hBkg->Draw("hist same");
	  hBkg->Fit("gaus");
	  auto c4 = new TCanvas("c4", "c4", 1000, 500);
	  c4->Divide(2,1);
	  c4->cd(1);
	  hAdcTime->Reset();   // We don't need this anymore
	  hAdcTime->SetAxisRange(0,1000,"Y");
	  hAdcTime->Draw();
	  for ( int colN=0; colN<32; colN++) {
	    for ( int rowN=0;rowN<8;rowN++) {
	      if ( label[colN][rowN] == 0 )  {
		
		if ( OddOrEven ( pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)) == 0)
		  hADCcollection[colN][rowN]->DrawCopy("same hist");
		
	      }
	    }
	  }
	  c4->cd(2);
	  hAdcTime->Draw();
          for ( int colN=0; colN<32; colN++) {
	    for ( int rowN=0;rowN<8;rowN++) {
              if ( label[colN][rowN] == 0 )  {
                if ( OddOrEven ( pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)) == 1)
                  hADCcollection[colN][rowN]->DrawCopy("same hist");

	      }
            }
          }
	  
	  
	  
	  auto c5 = new TCanvas("c5", "c5", 1000, 500);
	  c5->Divide(2,1);
	  c5->cd(1);
	  
	  for ( int colN=0; colN<32; colN++) {
	    for ( int rowN=0;rowN<8;rowN++) {
	      if ( label[colN][rowN] == 0 )  {
		
		if ( OddOrEven ( pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)) == 0) {
		  TH1F* htemp = (TH1F*)hADCcollection[colN][rowN]->Clone("htemp");
		  htemp->Scale( 100000. / htemp->Integral());
		  for ( int timeBuc = 0 ; timeBuc< 512; timeBuc++) {
		    BkgProfileEven->Fill ( timeBuc+1, htemp->GetBinContent( timeBuc+1) );
		  }
		  delete htemp;
		}
		
	      }
	    }
	  }
	  BkgProfileEven->Draw("colz");

	  c5->cd(2);
	  for ( int colN=0; colN<32; colN++) {
	    for ( int rowN=0;rowN<8;rowN++) {
	      if ( label[colN][rowN] == 0 )  {
		
		if ( OddOrEven ( pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)) == 1) {
		  TH1F* htemp = (TH1F*)hADCcollection[colN][rowN]->Clone("htemp");
		  htemp->Scale( 100000./ htemp->Integral());
		  for ( int timeBuc = 0 ; timeBuc< 512; timeBuc++) {
		    BkgProfileOdd->Fill ( timeBuc+1, htemp->GetBinContent( timeBuc+1) );
		  }	   
		}
		
	      }
	    }
	  }
	  BkgProfileOdd->Draw("colz");
	  

	// Remove the 5 max and 5 min ;
	vector<int> even_vec;
	vector<int> odd_vec;
	for ( int timeBuc = 0 ; timeBuc< 512; timeBuc++) {

	  even_vec.clear();
	  odd_vec.clear();
	  for ( int adcBin = 1;  adcBin <=nRelAdcBins ; adcBin++) { 
	    for ( int jj = 0 ; jj < BkgProfileEven->GetBinContent ( timeBuc+1, adcBin ) ;  jj++) {
	      even_vec.push_back( adcBin ) ;
	    }
	    
	    for ( int jj = 0 ; jj < BkgProfileOdd->GetBinContent ( timeBuc+1, adcBin ) ;  jj++) {
	      odd_vec.push_back( adcBin ) ;
	    }
	    
	  }

	  // Sort out
	  if ( even_vec.size() != 0 ) {
	    std::sort( even_vec.begin(), even_vec.end());
	    for ( int ii=3 ; ii< even_vec.size() -3; ii++) {
	      BkgProfileEvenCorr->Fill (  timeBuc+1, even_vec.at(ii)); 
	    }
	  }
	  
	  if ( odd_vec.size() != 0 ) {
	    std::sort( odd_vec.begin(), odd_vec.end() );
	    //	    std::sort( odd_vec.begin(), odd_vec.begin() + odd_vec.size());
	    for ( int ii=3 ; ii< odd_vec.size() -3; ii++) {
	      BkgProfileOddCorr->Fill (  timeBuc+1, odd_vec.at(ii)); 
	    }
	  }
	}
	
	TCanvas* c6 = new TCanvas("c6", "c6", 1000, 500);
	c6->Divide(2,1);
	c6->cd(1);
	BkgProfileEvenCorr->Draw("colz") ;
	c6->cd(2);
	BkgProfileOddCorr->Draw("colz") ;

	TCanvas* c7 = new TCanvas("c7", "c7", 1000, 500);
	c7->Divide(2,1);
	c7->cd(1);
	TH1F* finalBkgEven = (TH1F*) BkgProfileEvenCorr->ProfileX()->ProjectionX();
	finalBkgEven->Draw("e");
	
	c7->cd(2);
	TH1F* finalBkgOdd = (TH1F*) BkgProfileOddCorr->ProfileX()->ProjectionX();
	finalBkgOdd->Draw("e");

	TCanvas* c8 = new TCanvas("c8", "c8", 500, 500);
	for ( int colN=0; colN<32; colN++) {
	  for ( int rowN=0;rowN<8;rowN++) {
	    TH1F* hbkg;
	    if ( OddOrEven ( pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)) == 0) 
	      hbkg = finalBkgEven;
	    else if  ( OddOrEven ( pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)) == 1)
	      hbkg = finalBkgOdd;
	    else
	      continue;
	    
	    hbkg->SetLineColor(2);
	    hbkg->SetMarkerColor(2);
	    hbkg->Scale( hADCcollection[colN][rowN]->Integral(400,500) / hbkg->Integral(400,500) );

	    hADCcollectionBkgSubt[colN][rowN] = (TH1F*) hADCcollection[colN][rowN]->Clone(Form("adcBkgSubt_col%d_row%d",colN,rowN));
	    hADCcollectionBkgSubt[colN][rowN]->Add( hbkg, -1);
	    hADCcollectionBkgSubt[colN][rowN]->SetLineColor(4);

	    hADCcollection[colN][rowN]->SetAxisRange(-200,2000,"Y");
	    hADCcollection[colN][rowN]->SetLineWidth(2);
	    hbkg->SetLineWidth(2);
	    hADCcollectionBkgSubt[colN][rowN]->SetLineWidth(2);
	    hADCcollection[colN][rowN]->Draw();
	    hbkg->Draw("same hist");
	    hADCcollectionBkgSubt[colN][rowN]->Draw("same hist");
	    trimFirstFive( hADCcollectionBkgSubt[colN][rowN] ) ;  // We should clean up the first 5 time bucket due to some bugs. 
	    //	    hADCcollectionBkgSubt[colN][rowN]->Fit("fitFcn");
	    //	    fitFcn->SetLineColor(7);
	    //	    fitFcn->SetLineStyle(7);
	    //	    fitFcn->Draw("same");
	    
	    TLine* t1 = new TLine(0,0,512,0);
	    t1->SetLineStyle(7);
	    t1->Draw();
	    
	    
	    //	    c8->SaveAs(Form("./figures/BkgSubt_Event%d_%d_%d.png", nEvents,  pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)));
	  }}


	hPolyFinal->ClearBinContents();

	nhitsForTree =0 ;
	for ( int colN=0; colN<32; colN++) {
          for ( int rowN=0;rowN<8;rowN++) {
	    float max = hADCcollectionBkgSubt[colN][rowN]->GetBinContent(  hADCcollectionBkgSubt[colN][rowN]->GetMaximumBin());
	    //	    cout << " max = " << max << endl;

	    // Fill the tree
	    if ( max >= threshold2) {
	      adc[nhitsForTree] = max; 
	      time[nhitsForTree] =  hADCcollectionBkgSubt[colN][rowN]->GetMaximumBin(); 
	      x[nhitsForTree] = rowN * 100. / 8.; 
	      y[nhitsForTree] = colN * 100. / 32.; 
	      hPolyFinal->Fill(pMap.GetX(colN), pMap.GetY(colN, rowN), max);
	      
	      nhitsForTree++;
	    }
	    else
	      max =1 ;
	    
	    
	    
	  }}
	treeOut->Fill(); // fill the tree
      
	
	
	TCanvas* c9 = new TCanvas("c9", "c9", 700, 700);
	hPolyFinal->SetStats(0);
	hPolyFinal->SetAxisRange(1, 5000, "Z");
	hPolyFinal->SetTitle(Form("Event #%d;Y [mm];Z [mm]", nEvents));
	hPolyFinal->Draw("colz");
	c9->SetRightMargin(0.12);
	c9->Update();
	gPad->SetLogz();
	if (nHits > 4) {
	  c9->SaveAs(Form("./figures/trackFinal_%d.png", nEvents));
	}
	
	

	
      }
	
	nEvents++;
    }
    
    
    inFile.close();


    TFile* fout = new TFile("treeOfHits.root","recreate");
    treeOut->Write();
    fout->Close();
    
    
}
//double BackgroundFcn(double *x, double *par) {
//return par[0];
//}
double SignalFcn(double *x, double *par) {
  if (x[0] > par[1]) {
    return par[0] * TMath::Power((x[0] - par[1]), 3) * TMath::Exp(-par[2] * (x[0] - par[1]));
  } else {
    return 0;
  }
}
// double TotalFcn(double *x, double *par) {
//     return BackgroundFcn(x, par) + SignalFcn(x, &par[1]);
// }
// void DynamicExec() {
//     TObject *select = gPad->GetSelected();
//     if (!select) return;
//     if (!select->InheritsFrom(TH2::Class())) return;
//     TH2Poly *h = (TH2Poly *)select;
//     Float_t x = gPad->AbsPixeltoX(gPad->GetEventX());
//     Float_t y = gPad->AbsPixeltoY(gPad->GetEventY());
//     static int colN_save, rowN_save;
//     int colN, rowN;
//     colN = (h->FindBin(x, y) - 1) / 16;
//     rowN = (h->FindBin(x, y) - 1) % 16;
//     if ((colN_save == colN) && (rowN_save == rowN)) return;
//     colN_save = colN;
//     rowN_save = rowN;
//     TVirtualPad *padsav = gPad;
//     TCanvas *cExec = (TCanvas *)gROOT->FindObject("cExec");
//     if (cExec) {
//         // cExec->GetPrimitive("hExec")->Delete();
//         gROOT->FindObject("hExec")->Delete();
//         gROOT->FindObject("fitExec")->Delete();
//     } else {
//         cExec = new TCanvas("cExec", "", 700, 700);
//     }
//     auto *hExec = new TH1S("hExec", "", 512, 0, 512);
//     auto fitExec = new TF1("fitExec", TotalFcn, 0, 512, 4);
//     int meanValue = 0;
//     for (int buckIdx = 0; buckIdx < 512; buckIdx++) {
//         meanValue += Adc[colN][rowN][buckIdx];
//         hExec->SetBinContent(buckIdx + 1, Adc[colN][rowN][buckIdx]);
//     }
//     meanValue = TMath::Nint(meanValue / 512.);
//     int maxValue = hExec->GetBinContent(hExec->GetMaximumBin());
//     int maxBin = hExec->GetMaximumBin();
//     if (maxValue < meanValue + 50) {
//         cExec->cd();
//         hExec->SetStats(0);
//         hExec->SetAxisRange(0, 4500, "Y");
//         hExec->SetTitle(Form("(column, row) = (%d, %d), Max value : %d, Mean Value : %d;Time bucket;ADC (ch)", colN, rowN, maxValue, meanValue));
//         hExec->Draw();
//         cExec->Update();
//         padsav->cd();
//         return;
//     }
//     if (maxValue < 1000) {
//         fitExec->SetParameters(450, 0.1, maxBin, 0.1);
//         fitExec->SetParLimits(1, 0.001, 0.5);
//     } else if (maxValue < 2000) {
//         `fitExec->SetParameters(450, 1, maxBin, 0.1);
//         fitExec->SetParLimits(1, 0.3, 1.5);
//     } else {
//         fitExec->SetParameters(450, 2, maxBin, 0.1);
//         fitExec->SetParLimits(1, 1, 3);
//     }
//     fitExec->SetParLimits(0, 200, 600);
//     fitExec->SetParLimits(2, maxBin - 30, maxBin);
//     fitExec->SetParLimits(3, 0, 0.12);
//     cExec->cd();
//     hExec->Fit("fitExec", "Q");
//     hExec->SetStats(0);
//     hExec->SetAxisRange(0, 4500, "Y");
//     hExec->SetTitle(Form("(column, row) = (%d, %d), Max value : %d, Mean Value : %d;Time bucket;ADC (ch)", colN, rowN, maxValue, meanValue));
//     hExec->Draw();
//     cExec->Update();
//     double par[4];
//     fitExec->GetParameters(par);
//     cout << Form("Fit function f(x) = %lf + %lf * (x - %lf)^3 * exp[-%lf * (x - %lf)]", par[0], par[1], par[2], par[3], par[2]) << endl;
//     cout << "Reduced Chi square : " << fitExec->GetChisquare() / fitExec->GetNDF() << endl;
//     cout << "Integral Value : " << 2 * par[1] / TMath::Power(par[3], 3) << endl;
//     padsav->cd();
// }
