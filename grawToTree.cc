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

double bkgFcn(double *x, double *par) {
    return par[0];
}
double sigFcn(double *x, double *par) {
    return x[0] >= par[2] ? par[0] * TMath::Power(x[0] - par[2], 3) * TMath::Exp(-par[1] * (x[0] - par[2])) : 0;
}
double fitFunc(double *x, double *par) {
    return bkgFcn(x, par) + sigFcn(x, &par[1]);
}
bool isEven(int agetIdx, int chanIdx);

void grawToTree() {

  bool isDebugMode = true ;   // Save all supplemental figures 
  int threshold1 = 500 ; //   If the max ADC is smaller than threshold1, we assume that channel is background channel

    
  std::ifstream fList("files.txt");
  
    auto hSignal = new TH1F("hSignal", "", 511, 1, 512);
    auto htemp = (TH1F*)hSignal->Clone("htemp");  // will be used for temporary histogrmas
    TH1F* hSignalArr[4][68];  // Array for hSignal for all agets and channels
    for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
      for (int chanIdx = 0; chanIdx < 68; chanIdx++) {
	hSignalArr[agetIdx][chanIdx] = (TH1F*)hSignal->Clone(Form("hSignal_agetId%d_chanId%d",agetIdx,chanIdx));
      }
    }

    TH1F* hBkgTemplateOdd = (TH1F*)hSignal->Clone("hBkgTemplateOdd");
    TH1F* hBkgTemplateEven = (TH1F*)hSignal->Clone("hBkgTemplateEven");
	    
    
    TH1F *hBgkOddChanAget[4];
    TH1F *hBgkEvenChanAget[4];
    for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
        hBgkOddChanAget[agetIdx] = new TH1F(Form("hBgkOddChanAget%d", agetIdx), "", 511, 1, 512);
        hBgkEvenChanAget[agetIdx] = new TH1F(Form("hBgkEvenChanAget%d", agetIdx), "", 511, 1, 512);
    }
    auto hPolyPad = new TH2Poly("hPolyPad", "", -3.5625, 100.4375, -8.25, 95.75);
    auto hCharge = new TH1F("hCharge", "GEM_V = 310V;Pulse height(Channel);Counts", 150, 300, 1800);
    auto fFit = new TF1("fFit", fitFunc, 1, 512, 4);
    auto cvsSignal = new TCanvas("cvsSignal", "", 700, 700);
    auto cvsPad = new TCanvas("cvsPad", "", 700, 700);

    auto cvsAllChan = new TCanvas("cvsAllChan", "", 800, 800);    // Histograms for all channels will be drawn here
    cvsAllChan->Divide(2,2);    // for 4 aget channels
    auto cvsBkgChan = new TCanvas("cvsBkgChan", "", 800, 800);    // Histograms for all background channels (after normalization) will be drawn here
    cvsBkgChan->Divide(2,3);    // for odd and even channels
    
    // Histograms for background monitoring:
    int nRelAdcBins= 300;
    TH2F *BkgProfileOdd = new TH2F("BkgProfileOdd",";time bucket; Normalized ADC",511,1,512,nRelAdcBins,0,nRelAdcBins); // 2d histogram for time x ADC 
    TH2F *BkgProfileEven = new TH2F("BkgProfileEven",";time bucket; Normalized ADC",511,1,512,nRelAdcBins,0,nRelAdcBins);
    TH2F *BkgProfileOddCorr = new TH2F("BkgProfileOddCorr",";time bucket; Normalized ADC",511,1,512,nRelAdcBins,0,nRelAdcBins);
    TH2F *BkgProfileEvenCorr = new TH2F("BkgProfileEvenCorr",";time bucket; Normalized ADC",511,1,512,nRelAdcBins,0,nRelAdcBins);
    
    
    DataFrame frame;
    PadMap pMap;
    pMap.BuildPad(hPolyPad);

    while (!fList.eof()) {
        std::string fName;
        std::getline(fList, fName);
        if (fList.eof()) break;
        frame.OpenGrawFile(fName.c_str());
        int eventIdx;
        while (frame.Decode()) {
            eventIdx = frame.GetEventIdx();
            if (eventIdx % 10 == 0) {
                std::cout << eventIdx << std::endl;
            }
            hPolyPad->ClearBinContents();

            BkgProfileOdd->Reset();
            BkgProfileEven->Reset();
            BkgProfileOddCorr->Reset();
	    BkgProfileEvenCorr->Reset();
		    
	    for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
	      for (int chanIdx = 0; chanIdx < 68; chanIdx++) {
		hSignal->Reset();  // Initiation! 
		if (frame.IsFPNChannel(chanIdx)) continue;
		for (int buckIdx = 1; buckIdx < 512; buckIdx++) {
		  hSignal->SetBinContent(buckIdx, frame.GetADC(agetIdx, chanIdx, buckIdx));
		}
		
		// Background estimation procedure: 
		// step0 : copy the hSignal into the hSignalarrays
		hSignalArr[agetIdx][chanIdx]->Reset();
		hSignalArr[agetIdx][chanIdx]->Add(hSignal);

		// step1 : Select only background channels 
		bool isBkgChan = false; 
		int maxBin = hSignalArr[agetIdx][chanIdx]->GetMaximumBin();
		int maxVal = hSignalArr[agetIdx][chanIdx]->GetBinContent(maxBin);
		if (  maxVal < threshold1 ) 
		  isBkgChan = true; 

		// step2 : Normalize each histogram to fit the region of timebuck 0 ~ 100 ) 
		if ( isBkgChan == true) {

		  htemp->Reset();
		  htemp->Add(hSignalArr[agetIdx][chanIdx]);
		  htemp->Scale( 5000./ htemp->Integral(1,50)) ;
		  
		  for (int buckIdx = 1; buckIdx < 512; buckIdx++) {
		    if ( isEven(agetIdx, chanIdx) ) 
		      BkgProfileEven->Fill ( buckIdx, htemp->GetBinContent(buckIdx) );
		    else  
		      BkgProfileOdd->Fill ( buckIdx, htemp->GetBinContent(buckIdx) );
		    }		    
		}


		

	      }
	    }
	    // Remove the 5 maxs and 5 mins fomr BkgProfileEven and BkgProfileOdd
	    vector<int> even_vec;
	    vector<int> odd_vec;
	    for (int buckIdx = 1; buckIdx < 512; buckIdx++) {
	      even_vec.clear();
	      odd_vec.clear();
	      for ( int adcBin = 1;  adcBin <=nRelAdcBins ; adcBin++) { 
		for ( int jj = 0 ; jj < BkgProfileEven->GetBinContent ( buckIdx, adcBin ) ;  jj++) {
		  even_vec.push_back( adcBin ) ;
		}
		
		for ( int jj = 0 ; jj < BkgProfileOdd->GetBinContent ( buckIdx, adcBin ) ;  jj++) {
		  odd_vec.push_back( adcBin ) ;
		}
		
	      }
	    
	      // Sort out 
	      if ( even_vec.size() != 0 ) {
		std::sort( even_vec.begin(), even_vec.end());
		for ( int ii=3 ; ii< even_vec.size() -3 ; ii++) {
		  BkgProfileEvenCorr->Fill ( buckIdx, even_vec.at(ii)); 
		}
	      }
	      else {
		cout << "Not enough background channels!"  << endl;
		return ;
	      }	      
	      
	      // Fill the corrected histogram
	      if ( odd_vec.size() != 0 ) {
		std::sort( odd_vec.begin(), odd_vec.end() );
		//	    std::sort( odd_vec.begin(), odd_vec.begin() + odd_vec.size());
		for ( int ii=3 ; ii< odd_vec.size() -3 ;  ii++) {
		  BkgProfileOddCorr->Fill (  buckIdx, odd_vec.at(ii)); 
		}
	      }
	      else {
		cout << "Not enough background channels!"  << endl;
		return ;
	      }	      
	      
	    }

	    hBkgTemplateOdd->Reset();
	    hBkgTemplateEven->Reset();
	    hBkgTemplateOdd->Add( (TH1F*)BkgProfileOddCorr->ProfileX()->ProjectionX());
	    hBkgTemplateEven->Add( (TH1F*)BkgProfileEvenCorr->ProfileX()->ProjectionX());
	    
		
	    
	    // Check if the hisgroams are well made
	    if ( isDebugMode) {
	      for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
		cvsAllChan->cd(agetIdx+1);
		htemp->Reset();
		htemp->SetAxisRange(0,2000,"Y");
		htemp->Draw();
		for (int chanIdx = 0; chanIdx < 68; chanIdx++) {
		  hSignalArr[agetIdx][chanIdx]->Draw("same");
		}
	      }
	      cvsAllChan->SaveAs(Form("./figureDebug/cvsAllChan_%05d.png", eventIdx)); 
	      
	      cvsBkgChan->cd(1);
	      BkgProfileEven->Draw("colz");
	      cvsBkgChan->cd(2);
	      BkgProfileOdd->Draw("colz");

	      cvsBkgChan->cd(3);
	      BkgProfileEvenCorr->Draw("colz");
	      cvsBkgChan->cd(4);
	      BkgProfileOddCorr->Draw("colz");

	      cvsBkgChan->cd(5);
	      hBkgTemplateOdd->SetAxisRange(0,200,"Y");
	      hBkgTemplateOdd->Draw();
	      cvsBkgChan->cd(6);
	      hBkgTemplateEven->SetAxisRange(0,200,"Y");
	      hBkgTemplateEven->Draw();

	      cvsBkgChan->SaveAs(Form("./figureDebug/cvsBkgChan_%05d.png", eventIdx));

		
	    }
	    
	    
	    /*if (isEven(agetIdx, chanIdx)) {
	      hBgkEvenChanAget[agetIdx]->Add(hSignal, 1.0 / 32.);
	      } else {
	      hBgkOddChanAget[agetIdx]->Add(hSignal, 1.0 / 32.);
	      }
	    */
	    
	    // Reset monitoring histograms before entering the channel loop.
	    
	    
	    for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
	      for (int chanIdx = 0; chanIdx < 68; chanIdx++) {
		if (frame.IsFPNChannel(chanIdx)) continue;
		for (int buckIdx = 1; buckIdx < 512; buckIdx++) {
		  hSignal->SetBinContent(buckIdx, frame.GetADC(agetIdx, chanIdx, buckIdx));
                    }
		bool isSignalCand = false;
                    cvsSignal->cd();
                    if (isEven(agetIdx, chanIdx)) {
                        hSignal->Add(hBgkEvenChanAget[agetIdx], -1.0);
                        int maxBin = hSignal->GetMaximumBin();
                        int maxVal = hSignal->GetBinContent(maxBin);
                        double meanVal = hSignal->Integral(1, 50) / 50.;
                        // hSignal->SetTitle(Form("%lf", meanVal));
                        hSignal->SetAxisRange(-100, 1000, "Y");
                        hSignal->Draw();
                        if (maxBin < 450 && maxBin > 50 && maxVal - meanVal > 20) {
			  fFit->SetParameter(0,2);   // Should set initial parameters for fair fits.  Otherwise, the previous channel result will bias it. 
			  fFit->SetParameter(1,1000);
			  fFit->SetParameter(2,1000);
			  fFit->SetParameter(3,10);

			  fFit->SetParLimits(0, meanVal - 10, meanVal + 10);
                            fFit->SetParLimits(1, 0.001, 0.6);
                            fFit->SetParLimits(2, 0.048, 0.052);
                            fFit->SetParLimits(3, 0, maxBin);
                            fFit->SetParameters(meanVal, maxVal * TMath::Power((0.05 * TMath::E()) / 3., 3) + 0.02, 0.05, maxBin - 50);
                            fFit->SetRange(0, maxBin + 50);
                            hSignal->Fit("fFit", "Q");
                            gStyle->SetOptFit(1111);
                            isSignalCand = true;
                        }
                    } else {
                        hSignal->Add(hBgkOddChanAget[agetIdx], -1.0);
                        int maxBin = hSignal->GetMaximumBin();
                        int maxVal = hSignal->GetBinContent(maxBin);
                        double meanVal = hSignal->Integral(1, 50) / 50.;
                        hSignal->SetTitle(Form("%lf", meanVal));
                        hSignal->SetAxisRange(-100, 1000, "Y");
                        hSignal->Draw();
                        if (maxBin < 450 && maxBin > 50 && maxVal - meanVal > 20) {
			  fFit->SetParameter(0,2);   // Should set initial parameters for fair fits.  Otherwise, the previous channel result will bias it. 
			  fFit->SetParameter(1,1000);
			  fFit->SetParameter(2,1000);
			  fFit->SetParameter(3,10);
                            fFit->SetParLimits(0, meanVal - 10, meanVal + 10);
                            fFit->SetParLimits(1, 0.001, 0.6);
                            fFit->SetParLimits(2, 0.048, 0.052);
                            fFit->SetParLimits(3, 0, maxBin);
                            fFit->SetParameters(meanVal, maxVal * TMath::Power((0.05 * TMath::E()) / 3., 3) + 0.02, 0.05, maxBin - 50);
                            fFit->SetRange(0, maxBin + 50);
                            hSignal->Fit("fFit", "RQ");
                            gStyle->SetOptFit(1111);
                            isSignalCand = true;
                        }
                    }
		    if (isSignalCand) {
                        double aa = fFit->GetParameter(1);
                        double bb = fFit->GetParameter(3);
                        double cc = fFit->GetParameter(2);
                        double realMaxVal = aa * TMath::Power(3. / (cc * TMath::E()), 3);
                        if (1) {
                            int xCoor = pMap.GetX(agetIdx, chanIdx);
                            int yCoor = pMap.GetY(agetIdx, chanIdx);
                            hPolyPad->Fill(xCoor, yCoor, realMaxVal);
                            // hSignal->SetTitle(Form("%lf, %lf, %lf", aa, cc, realMaxVal));
                        }
			
			if ( isDebugMode) { 
			  hSignal->Draw();
			  cvsSignal->Update();
			  //			  cvsSignal->SaveAs(Form("./fitResults/signal_%05d_%d_%d.png", eventIdx, agetIdx, chanIdx));
			}
                    }
		}
            }
            cvsPad->cd();
            hPolyPad->SetAxisRange(0, 4000, "Z");
            hPolyPad->Draw("colz");
            hPolyPad->SetStats(0);
            cvsPad->Update();
            cvsPad->SetLogz();
            cvsPad->SaveAs(Form("./track/event%05d.png", eventIdx));
        }
        frame.CloseGrawFile();
    }
}
bool isEven(int agetIdx, int chanIdx) {
    if (agetIdx == 0) {  // 순서 바뀜
        if ((chanIdx < 11 || (chanIdx > 22 && chanIdx < 45) || chanIdx > 56) && chanIdx % 2 == 1) {
            return true;
        } else if (((chanIdx > 11 && chanIdx < 22) || (chanIdx > 45 && chanIdx < 56)) && chanIdx % 2 == 0) {
            return true;
        } else {
            return false;
        }
    } else {
        if ((chanIdx < 11 || (chanIdx > 22 && chanIdx < 45) || chanIdx > 56) && chanIdx % 2 == 0) {
            return true;
        } else if (((chanIdx > 11 && chanIdx < 22) || (chanIdx > 45 && chanIdx < 56)) && chanIdx % 2 == 1) {
            return true;
        } else {
            return false;
        }
    }
}
