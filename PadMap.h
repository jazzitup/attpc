#ifndef __PadMap_h
#define __PadMap_h

#include <cmath>
#include <map>
#include <vector>

#include "TH2Poly.h"

class PadMap {
private:
    int padID_;
    std::map<std::pair<int, int>, std::vector<int>> mapToID_;
    std::map<std::pair<int, int>, std::vector<int>> mapToPosition_;

public:
    PadMap(int padID = 1);
    void TransformChannelToPad(int &agetID, int &channelID, int &colN, int &rowN);
    // HoneyComb Pad: 0, Rectangle Pad: 1
    int GetAgetIdx(int colN, int rowN);
    int GetChanIdx(int colN, int rowN);
    // int GetFPNChannelID(int channelID);
    int GetColN(int agetID, int channelID);
    int GetRowN(int agetID, int channelID);
    float GetX(int colN);
    float GetY(int colN, int rowN);
    void BuildPad(TH2Poly *poly);
};

#endif