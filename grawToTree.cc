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
    std::ifstream fList("files.txt");

    auto hSignal = new TH1F("hSignal", "", 511, 1, 512);
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
            for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
                hBgkOddChanAget[agetIdx]->Reset("ICESM");
                hBgkEvenChanAget[agetIdx]->Reset("ICESM");
		hSignal->Reset();  // Initiation! 
		for (int chanIdx = 0; chanIdx < 68; chanIdx++) {
                    if (frame.IsFPNChannel(chanIdx)) continue;
                    for (int buckIdx = 1; buckIdx < 512; buckIdx++) {
                        hSignal->SetBinContent(buckIdx, frame.GetADC(agetIdx, chanIdx, buckIdx));
                    }
                    if (isEven(agetIdx, chanIdx)) {
                        hBgkEvenChanAget[agetIdx]->Add(hSignal, 1.0 / 32.);
                    } else {
                        hBgkOddChanAget[agetIdx]->Add(hSignal, 1.0 / 32.);
                    }
                }
            }
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
                        int minBin = hSignal->GetMinimumBin();
                        int maxBin = hSignal->GetMaximumBin();
                        int minVal = hSignal->GetBinContent(minBin);
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
                            hSignal->Fit("fFit", "");
                            gStyle->SetOptFit(1111);
                            isSignalCand = true;
                        }
                    } else {
                        hSignal->Add(hBgkOddChanAget[agetIdx], -1.0);
                        int minBin = hSignal->GetMinimumBin();
                        int maxBin = hSignal->GetMaximumBin();
                        int minVal = hSignal->GetBinContent(minBin);
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

			hSignal->Draw();
			cvsSignal->Update();
			cvsSignal->SaveAs(Form("./fitResults/signal_%05d_%d_%d.png", eventIdx, agetIdx, chanIdx));
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
