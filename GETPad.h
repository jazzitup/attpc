#ifndef __GETPad_h
#define __GETPad_h

#include <iostream>

#include "TH2Poly.h"

class GETPad {
private:
    const double dx_ = 2.625;
    const double dy_ = 12.;
    const double gap_ = 0.5;
    int xId_[4][68];
    int yId_[4][68];
    int agetId_[32][8];
    int chanId_[32][8];

public:
    GETPad();
    int GetAgetId(int xId, int yId);
    int GetChanId(int xId, int yId);
    int GetXId(int agetId, int chanId);
    int GetYId(int agetId, int chanId);
    double GetX(int agetId, int chanId);
    double GetY(int agetId, int chanId);
    void AddBins(TH2Poly *poly);
    bool IsDeadPad(int xId, int yId);
    bool IsEvenChan(int agetId, int chanId);
    bool IsDeadChan(int agetId, int chanId);
};
#endif
