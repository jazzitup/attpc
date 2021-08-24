#ifndef __GETAnalyzer_h
#define __GETAnalyzer_h

#include <algorithm>
#include <numeric>
#include <cmath>

#include "GETDecoder.h"

class GETAnalyzer {
   private:
    GETDecoder *decoder_;
    bool isFakeEvent_;
    bool isSparkEvent_;
    int nFakeEvent_;
    int nSparkEvent_;
    const int FPN_chanId_[4] = {11, 22, 45 , 56};
    double ADC_ANC_[4][68][512];
    double ADC_FPNC_[4][68][512];
    double noiseEvn_[4][512];
    double noiseOdd_[4][512];

   public:
    GETAnalyzer();
    void LinkToDecoder(GETDecoder *decoder);
    void RunEventChecker();
    void RunNoiseCanceller();
    double GetADC_ANC(int agetId, int chanId, int buckId);
    double GetADC_FPNC(int agetId, int chanId, int buckId);
    double GetNoise(int agetId, int buckId, bool useEvn);
    bool IsFakeEvent();
    bool IsSparkEvent();
    int GetNumberOfFakeEventBefore();
    int GetNumberOfSparkEventBefore();
    bool IsFPNChan(int chanId);
    bool IsEvenChan(int agetId, int chanId);
    bool IsDeadChan(int agetId, int chanId);
};
//
bool compareFirstByDescending(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return (a.first > b.first);
}
bool compareFirstByAscending(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return (a.first < b.first);
}
bool compareSecondByDescending(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return (a.second > b.second);
}
bool compareSecondByAscending(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return (a.second < b.second);
}
#endif