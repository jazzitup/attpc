#ifndef __PadMap_h
#define __PadMap_h

#include <cmath>
#include <map>
#include <vector>

#include "TH2Poly.h"

class PadMap {
private:
    std::map<std::pair<int, int>, std::pair<int, int>> padToAgetMap_;
    std::map<std::pair<int, int>, std::pair<int, int>> agetToPadMap_;

public:
    PadMap();
    void TransformChannelToPad(int &agetIdx, int &chanIdx, int &colIdx, int &rowIdx);
    int GetAgetIdx(int colIdx, int rowIdx);
    int GetChanIdx(int colIdx, int rowIdx);
    // int GetFPNChannelID(int channelID);
    int GetColIdx(int agetIdx, int chanIdx);
    int GetRowIdx(int agetIdx, int chanIdx);
    float GetX(int agetIdx, int chanIdx);
    float GetY(int agetIdx, int chanIdx);
    void BuildPad(TH2Poly *poly);
};

#endif