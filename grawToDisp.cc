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
// double SignalFcn(double *x, double *par);
// double TotalFcn(double *x, double *par);


using std::cout;
using std::endl;
void grawToDisp(std::string strInFileName = "RawFiles/cosmic_202104/CoBo_2021-05-13T19h00m11.474s_0000.graw") {
    std::ifstream inFile;
    inFile.open(strInFileName.c_str(), std::ios::in | std::ios::binary);
    if (!inFile.is_open()) {
        cout << "File open error : " << strInFileName << endl;
        exit(-1);
    }
    auto hPoly = new TH2Poly("Read-out Pad", "", -3.5625, 100.4375, -8.25, 95.75);
    auto c1 = new TCanvas("c1", "c1", 700, 700);
    auto c2 = new TCanvas("c2", "c2", 700, 700);
    auto *hAdcTime = new TH1F("hAdcTime", "", 512, 0, 512);
    // auto fitFcn = new TF1("fitFcn", TotalFcn, 0, 512, 4);
    DataFrame frame;
    PadMap pMap(1);
    pMap.BuildPad(hPoly);
    int nEvents = 0;
    while (!inFile.eof()) {
        frame.ReadHeader(inFile);
        if (inFile.eof()) {
            cout << "[EOF]" << endl;
            break;
        }
        frame.ReadItem(inFile);
        hPoly->ClearBinContents();
        int nHits = 0;
        cout << nEvents << endl;
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
                c1->cd();
                hAdcTime->SetStats(0);
                hAdcTime->SetAxisRange(0, 4500, "Y");
                c1->Update();
                hAdcTime->Draw();
                c1->SaveAs(Form("./figures/Event%d_%d_%d.png", nEvents, pMap.GetAgetIdx(colN, rowN), pMap.GetChanIdx(colN, rowN)));
                if (maxValue > 600) {
                    hPoly->Fill(pMap.GetX(colN), pMap.GetY(colN, rowN), maxValue);
                    nHits++;
                } else {
                    hPoly->Fill(pMap.GetX(colN), pMap.GetY(colN, rowN), 1);
                }
            }
        }
        c2->cd();
        c2->SetLogz();
        hPoly->SetStats(0);
        hPoly->SetAxisRange(1, 5000, "Z");
        hPoly->SetTitle(Form("Event #%d;Y [mm];Z [mm]", nEvents));
        hPoly->Draw("colz");
        c2->Update();
        if (nHits > 4) {
            c2->SaveAs(Form("./figures/track_%d.png", nEvents));
        }
        nEvents++;
    }
    inFile.close();
}
// double BackgroundFcn(double *x, double *par) {
//     return par[0];
// }
// double SignalFcn(double *x, double *par) {
//     if (x[0] > par[1]) {
//         return par[0] * TMath::Power((x[0] - par[1]), 3) * TMath::Exp(-par[2] * (x[0] - par[1]));
//     } else {
//         return 0;
//     }
// }
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
//         fitExec->SetParameters(450, 1, maxBin, 0.1);
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
