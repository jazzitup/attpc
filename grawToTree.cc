#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "GETAnalyzer.h"
#include "GETDecoder.h"
#include "GETPad.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "TPad.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

const double par3 = 2.685;
const double par4 = 16.01;
const double par1ToADC = TMath::Power(par3 * par4 / TMath::E(), par3);

double sigFcn(double *x, double *par) {
    return x[0] < par[2] ? par[0] : par[1] * TMath::Power(x[0] - par[2], par3) * TMath::Exp(-(x[0] - par[2]) / par4);
}
double fitFcn(double *x, double *par) {
    return sigFcn(x, par) < par[3] ? sigFcn(x, par) : par[3];
}
void grawToTree(int runId = 1) {
  
    bool isDebugMode = 0 ;   // Save the performance plots
    TFile *fout = new TFile(Form("treeOfHits_run%d_v6_2021Aug24.root", runId), "recreate");
    GETDecoder decoder;
    if ( runId == 0 ) { // test run
      decoder.OpenFromList("fileList/files_muon_run1.txt");
    }
    else if ( runId <=9) {
      decoder.OpenFromList(Form("fileList/files_muon_run%d.txt", runId));
    }
    else {
      cout << " No event for runId = " << runId << endl;
      return;
    }
    
    GETAnalyzer analyzer;
    analyzer.LinkToDecoder(&decoder);
    GETPad pad;
    //
    auto treeOut = new TTree("treeOfHits", "a Tree with hits");  // output tree
    int eventId;
    bool isSpark;
    double eventTime, diffTime;  // (10 ns unit)
    int nHits;
    int runNumber = runId;
    double xTree[256];
    double yTree[256];
    int xIdTree[256];
    int yIdTree[256];
    int agetIdTree[256];
    int chanIdTree[256];
    double timeTree[256];
    double adcTree[256];
    treeOut->Branch("run", &runNumber, "run/I");
    treeOut->Branch("eventId", &eventId, "eventId/I");
    treeOut->Branch("isSpark", &isSpark, "isSpark/O");
    treeOut->Branch("eventTime", &eventTime, "eventTime/D");
    treeOut->Branch("diffTime", &diffTime, "diffTime/D");
    treeOut->Branch("nHits", &nHits, "nHits/I");
    treeOut->Branch("x", xTree, "x[nHits]/D");
    treeOut->Branch("y", yTree, "y[nHits]/D");
    treeOut->Branch("xId", xIdTree, "xId[nHits]/I");
    treeOut->Branch("yId", yIdTree, "yId[nHits]/I");
    treeOut->Branch("agetId", agetIdTree, "agetId[nHits]/I");
    treeOut->Branch("chanId", chanIdTree, "chanId[nHits]/I");
    treeOut->Branch("time", timeTree, "time[nHits]/D");
    treeOut->Branch("adc", adcTree, "adc[nHits]/D");
    //
    auto hTemplate = new TH1D("hTemplate", "", 500, 0.5, 500.5);
    auto fFitFcn = new TF1("fFitFcn", fitFcn, 1, 500, 4);
    auto c1 = new TCanvas("c1", "", 2000, 1400);
    double startTime, prevTime;
    //    int counter_for_test = 0;  
    while (decoder.Run()) {
      //      counter_for_test++;
      analyzer.RunEventChecker();
        if (analyzer.IsFakeEvent()) continue;
        eventId = decoder.GetEventId() - analyzer.GetNumberOfFakeEventBefore();
        if (eventId == 0) {
            startTime = (double)decoder.GetEventTime() / 1.E8;
            eventTime = 0;
            diffTime = 0;
            prevTime = 0;
        } else {
            eventTime = (double)decoder.GetEventTime() / 1.E8 - startTime;
            diffTime = eventTime - prevTime;
            prevTime = eventTime;
        }
        std::cout << eventId << std::endl;
        nHits = 0;
        if (analyzer.IsSparkEvent()) {
            isSpark = true;
            treeOut->Fill();
            continue;
        }
        isSpark = false;
        analyzer.RunNoiseCanceller();
        for (int agetId = 0; agetId < 4; agetId++) {
            for (int chanId = 0; chanId < 68; chanId++) {
                if (analyzer.IsFPNChan(chanId) || analyzer.IsDeadChan(agetId, chanId)) continue;
                auto hSignal = (TH1D *)hTemplate->Clone("hSignal");
                for (int buckId = 1; buckId <= 500; buckId++) {
                    hSignal->SetBinContent(buckId, analyzer.GetADC_ANC(agetId, chanId, buckId));
                }
                int maxBin = hSignal->GetMaximumBin();
                double maxADC = hSignal->GetBinContent(maxBin);
                if (isDebugMode) {
		  c1->cd();
		  hSignal->SetAxisRange(-100, 4000, "Y");
		  hSignal->SetLineColor(kBlue);
		  hSignal->Draw("hist");
		}
		//		if(maxADC < 3400 && maxADC > 50) {
		if (maxADC > 50) {
		  if (maxADC > 3400) {
		    fFitFcn->SetParameters(0, maxADC / par1ToADC, maxBin - par3 * par4, maxADC + 100);
		    fFitFcn->SetParLimits(1, (maxADC - 5) / par1ToADC, (maxADC + 5000) / par1ToADC);
		    fFitFcn->SetParLimits(2, maxBin - 100, maxBin);
		  } else {
		    fFitFcn->SetParameters(0, maxADC / par1ToADC, maxBin - par3 * par4, maxADC + 100);
		    fFitFcn->SetParLimits(1, (maxADC - 5) / par1ToADC, (maxADC + 5) / par1ToADC);
		    fFitFcn->SetParLimits(2, maxBin - par3 * par4 - 10, maxBin - par3 * par4 + 10);
		  }
		  fFitFcn->SetParLimits(0, -10, 10);
		  fFitFcn->SetParLimits(3, maxADC + 10, maxADC + 100);
		  fFitFcn->SetLineColor(kRed);
		  fFitFcn->SetLineWidth(2);
		  fFitFcn->SetRange(maxBin - 25, maxBin + 50);
		  hSignal->Fit(fFitFcn, "QNR");
		  double pars[4];
		  fFitFcn->GetParameters(pars);
		    if (isDebugMode) {
		      hSignal->SetTitle(Form("%.1lf + %.3lf x (x[0] - %.1lf)^2.685 x exp(-(x[0] - %.1lf) / 16.01);TimeBucket (10 ns);ADC", pars[0], pars[1], pars[2], pars[2]));
		      fFitFcn->SetRange(1, 500);
		      fFitFcn->Draw("same");
		      c1->Update();
		      c1->SaveAs(Form("./signal/event_run%02d_%05d_%d_%d.png", runId, eventId, agetId, chanId));
		    }
		    agetIdTree[nHits] = agetId;
		    chanIdTree[nHits] = chanId;
                    xIdTree[nHits] = pad.GetXId(agetId, chanId);
                    yIdTree[nHits] = pad.GetYId(agetId, chanId);
                    xTree[nHits] = pad.GetX(agetId, chanId);
                    yTree[nHits] = pad.GetY(agetId, chanId);
                    timeTree[nHits] = pars[2];
                    adcTree[nHits] = pars[1] * par1ToADC;
                    nHits++;
		}
                hSignal->Delete();
	    }
        }
        treeOut->Fill();
    }
    treeOut->Write();
    fout->Close();
}
