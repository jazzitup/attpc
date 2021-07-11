#include "PadMap.h"
//
PadMap::PadMap() {
    std::pair<int, int> agetPair;
    std::pair<int, int> padPair;
    int agetIdx, chanIdx, colIdx, rowIdx;
    for (agetIdx = 0; agetIdx < 4; agetIdx++) {
        for (chanIdx = 0; chanIdx < 68; chanIdx++) {
            if (chanIdx == 11 || chanIdx == 22 || chanIdx == 45 || chanIdx == 56) continue;  // FPN Channel
            TransformChannelToPad(agetIdx, chanIdx, colIdx, rowIdx);
            agetPair = std::make_pair(agetIdx, chanIdx);
            padPair = std::make_pair(colIdx, rowIdx);
            this->padToAgetMap_[padPair] = agetPair;
            this->agetToPadMap_[agetPair] = padPair;
        }
    }
}
int PadMap::GetAgetIdx(int colIdx, int rowIdx) {
    return padToAgetMap_[std::make_pair(colIdx, rowIdx)].first;
}
int PadMap::GetChanIdx(int colIdx, int rowIdx) {
    return padToAgetMap_[std::make_pair(colIdx, rowIdx)].second;
}
int PadMap::GetColIdx(int agetIdx, int chanIdx) {
    return agetToPadMap_[std::make_pair(agetIdx, chanIdx)].first;
}
int PadMap::GetRowIdx(int agetIdx, int chanIdx) {
    return agetToPadMap_[std::make_pair(agetIdx, chanIdx)].second;
}
float PadMap::GetX(int agetIdx, int chanIdx) {
    return 3.125 * agetToPadMap_[std::make_pair(agetIdx, chanIdx)].first;
}
float PadMap::GetY(int agetIdx, int chanIdx) {
    return 12.5 * agetToPadMap_[std::make_pair(agetIdx, chanIdx)].second;
}
int PadMap::GetXId(int agetIdx, int chanIdx) {
    return agetToPadMap_[std::make_pair(agetIdx, chanIdx)].first;
}
int PadMap::GetYId(int agetIdx, int chanIdx) {
    return agetToPadMap_[std::make_pair(agetIdx, chanIdx)].second;
}



void PadMap::BuildPad(TH2Poly *poly) {
    // Add the bins
    double x[4], y[4];
    // Set Rectangle pattern parameters
    double par[3];
    par[0] = 2.625;
    par[1] = 12;
    par[2] = 0.5;  // x-spacing, y-spacing, gap spacing
    double xCenter = 0, yCenter = 0;
    for (int nColumn = 0; nColumn < 32; nColumn++) {
        for (int nRow = 0; nRow < 8; nRow++) {
            x[0] = xCenter - 0.5 * par[0];
            y[0] = yCenter - 0.5 * par[1];
            x[1] = x[0];
            y[1] = yCenter + 0.5 * par[1];
            x[2] = xCenter + 0.5 * par[0];
            y[2] = y[1];
            x[3] = x[2];
            y[3] = y[0];

            poly->AddBin(4, x, y);

            // Go up
            yCenter += par[1] + par[2];
        }
        xCenter += par[0] + par[2];
        yCenter = 0;
    }
}
void PadMap::TransformChannelToPad(int &agetIdx, int &chanIdx, int &colIdx, int &rowIdx) {
    if (agetIdx == 0) {
        switch (chanIdx) {
            case 0:
                colIdx = 16;
                rowIdx = 3;
                break;
            case 1:
                colIdx = 17;
                rowIdx = 3;
                break;
            case 2:
                colIdx = 18;
                rowIdx = 3;
                break;
            case 3:
                colIdx = 19;
                rowIdx = 3;
                break;
            case 4:
                colIdx = 20;
                rowIdx = 3;
                break;
            case 5:
                colIdx = 21;
                rowIdx = 3;
                break;
            case 6:
                colIdx = 22;
                rowIdx = 3;
                break;
            case 7:
                colIdx = 23;
                rowIdx = 3;
                break;
            case 8:
                colIdx = 24;
                rowIdx = 3;
                break;
            case 9:
                colIdx = 25;
                rowIdx = 3;
                break;
            case 10:
                colIdx = 26;
                rowIdx = 3;
                break;
            case 12:
                colIdx = 27;
                rowIdx = 3;
                break;
            case 13:
                colIdx = 28;
                rowIdx = 3;
                break;
            case 14:
                colIdx = 29;
                rowIdx = 3;
                break;
            case 15:
                colIdx = 30;
                rowIdx = 3;
                break;
            case 16:
                colIdx = 31;
                rowIdx = 3;
                break;
            case 17:
                colIdx = 16;
                rowIdx = 2;
                break;
            case 18:
                colIdx = 17;
                rowIdx = 2;
                break;
            case 19:
                colIdx = 18;
                rowIdx = 2;
                break;
            case 20:
                colIdx = 19;
                rowIdx = 2;
                break;
            case 21:
                colIdx = 20;
                rowIdx = 2;
                break;
            case 23:
                colIdx = 21;
                rowIdx = 2;
                break;
            case 24:
                colIdx = 22;
                rowIdx = 2;
                break;
            case 25:
                colIdx = 23;
                rowIdx = 2;
                break;
            case 26:
                colIdx = 24;
                rowIdx = 2;
                break;
            case 27:
                colIdx = 25;
                rowIdx = 2;
                break;
            case 28:
                colIdx = 26;
                rowIdx = 2;
                break;
            case 29:
                colIdx = 27;
                rowIdx = 2;
                break;
            case 30:
                colIdx = 28;
                rowIdx = 2;
                break;
            case 31:
                colIdx = 29;
                rowIdx = 2;
                break;
            case 32:
                colIdx = 30;
                rowIdx = 2;
                break;
            case 33:
                colIdx = 31;
                rowIdx = 2;
                break;
            case 34:
                colIdx = 16;
                rowIdx = 1;
                break;
            case 35:
                colIdx = 17;
                rowIdx = 1;
                break;
            case 36:
                colIdx = 18;
                rowIdx = 1;
                break;
            case 37:
                colIdx = 19;
                rowIdx = 1;
                break;
            case 38:
                colIdx = 20;
                rowIdx = 1;
                break;
            case 39:
                colIdx = 21;
                rowIdx = 1;
                break;
            case 40:
                colIdx = 22;
                rowIdx = 1;
                break;
            case 41:
                colIdx = 23;
                rowIdx = 1;
                break;
            case 42:
                colIdx = 24;
                rowIdx = 1;
                break;
            case 43:
                colIdx = 25;
                rowIdx = 1;
                break;
            case 44:
                colIdx = 26;
                rowIdx = 1;
                break;
            case 46:
                colIdx = 27;
                rowIdx = 1;
                break;
            case 47:
                colIdx = 28;
                rowIdx = 1;
                break;
            case 48:
                colIdx = 29;
                rowIdx = 1;
                break;
            case 49:
                colIdx = 30;
                rowIdx = 1;
                break;
            case 50:
                colIdx = 31;
                rowIdx = 1;
                break;
            case 51:
                colIdx = 16;
                rowIdx = 0;
                break;
            case 52:
                colIdx = 17;
                rowIdx = 0;
                break;
            case 53:
                colIdx = 18;
                rowIdx = 0;
                break;
            case 54:
                colIdx = 19;
                rowIdx = 0;
                break;
            case 55:
                colIdx = 20;
                rowIdx = 0;
                break;
            case 57:
                colIdx = 21;
                rowIdx = 0;
                break;
            case 58:
                colIdx = 22;
                rowIdx = 0;
                break;
            case 59:
                colIdx = 23;
                rowIdx = 0;
                break;
            case 60:
                colIdx = 24;
                rowIdx = 0;
                break;
            case 61:
                colIdx = 25;
                rowIdx = 0;
                break;
            case 62:
                colIdx = 26;
                rowIdx = 0;
                break;
            case 63:
                colIdx = 27;
                rowIdx = 0;
                break;
            case 64:
                colIdx = 28;
                rowIdx = 0;
                break;
            case 65:
                colIdx = 29;
                rowIdx = 0;
                break;
            case 66:
                colIdx = 30;
                rowIdx = 0;
                break;
            case 67:
                colIdx = 31;
                rowIdx = 0;
                break;

            default:
                colIdx = -99;
                rowIdx = -99;
                return;
        }
    } else if (agetIdx == 1) {
        switch (chanIdx) {
            case 0:
                colIdx = 16;
                rowIdx = 7;
                break;
            case 1:
                colIdx = 17;
                rowIdx = 7;
                break;
            case 2:
                colIdx = 18;
                rowIdx = 7;
                break;
            case 3:
                colIdx = 19;
                rowIdx = 7;
                break;
            case 4:
                colIdx = 20;
                rowIdx = 7;
                break;
            case 5:
                colIdx = 21;
                rowIdx = 7;
                break;
            case 6:
                colIdx = 22;
                rowIdx = 7;
                break;
            case 7:
                colIdx = 23;
                rowIdx = 7;
                break;
            case 8:
                colIdx = 24;
                rowIdx = 7;
                break;
            case 9:
                colIdx = 25;
                rowIdx = 7;
                break;
            case 10:
                colIdx = 26;
                rowIdx = 7;
                break;
            case 12:
                colIdx = 27;
                rowIdx = 7;
                break;
            case 13:
                colIdx = 28;
                rowIdx = 7;
                break;
            case 14:
                colIdx = 29;
                rowIdx = 7;
                break;
            case 15:
                colIdx = 30;
                rowIdx = 7;
                break;
            case 16:
                colIdx = 31;
                rowIdx = 7;
                break;
            case 17:
                colIdx = 16;
                rowIdx = 6;
                break;
            case 18:
                colIdx = 17;
                rowIdx = 6;
                break;
            case 19:
                colIdx = 18;
                rowIdx = 6;
                break;
            case 20:
                colIdx = 19;
                rowIdx = 6;
                break;
            case 21:
                colIdx = 20;
                rowIdx = 6;
                break;
            case 23:
                colIdx = 21;
                rowIdx = 6;
                break;
            case 24:
                colIdx = 22;
                rowIdx = 6;
                break;
            case 25:
                colIdx = 23;
                rowIdx = 6;
                break;
            case 26:
                colIdx = 24;
                rowIdx = 6;
                break;
            case 27:
                colIdx = 25;
                rowIdx = 6;
                break;
            case 28:
                colIdx = 26;
                rowIdx = 6;
                break;
            case 29:
                colIdx = 27;
                rowIdx = 6;
                break;
            case 30:
                colIdx = 28;
                rowIdx = 6;
                break;
            case 31:
                colIdx = 29;
                rowIdx = 6;
                break;
            case 32:
                colIdx = 30;
                rowIdx = 6;
                break;
            case 33:
                colIdx = 31;
                rowIdx = 6;
                break;
            case 34:
                colIdx = 16;
                rowIdx = 5;
                break;
            case 35:
                colIdx = 17;
                rowIdx = 5;
                break;
            case 36:
                colIdx = 18;
                rowIdx = 5;
                break;
            case 37:
                colIdx = 19;
                rowIdx = 5;
                break;
            case 38:
                colIdx = 20;
                rowIdx = 5;
                break;
            case 39:
                colIdx = 21;
                rowIdx = 5;
                break;
            case 40:
                colIdx = 22;
                rowIdx = 5;
                break;
            case 41:
                colIdx = 23;
                rowIdx = 5;
                break;
            case 42:
                colIdx = 24;
                rowIdx = 5;
                break;
            case 43:
                colIdx = 25;
                rowIdx = 5;
                break;
            case 44:
                colIdx = 26;
                rowIdx = 5;
                break;
            case 46:
                colIdx = 27;
                rowIdx = 5;
                break;
            case 47:
                colIdx = 28;
                rowIdx = 5;
                break;
            case 48:
                colIdx = 29;
                rowIdx = 5;
                break;
            case 49:
                colIdx = 30;
                rowIdx = 5;
                break;
            case 50:
                colIdx = 31;
                rowIdx = 5;
                break;
            case 51:
                colIdx = 16;
                rowIdx = 4;
                break;
            case 52:
                colIdx = 17;
                rowIdx = 4;
                break;
            case 53:
                colIdx = 18;
                rowIdx = 4;
                break;
            case 54:
                colIdx = 19;
                rowIdx = 4;
                break;
            case 55:
                colIdx = 20;
                rowIdx = 4;
                break;
            case 57:
                colIdx = 21;
                rowIdx = 4;
                break;
            case 58:
                colIdx = 22;
                rowIdx = 4;
                break;
            case 59:
                colIdx = 23;
                rowIdx = 4;
                break;
            case 60:
                colIdx = 24;
                rowIdx = 4;
                break;
            case 61:
                colIdx = 25;
                rowIdx = 4;
                break;
            case 62:
                colIdx = 26;
                rowIdx = 4;
                break;
            case 63:
                colIdx = 27;
                rowIdx = 4;
                break;
            case 64:
                colIdx = 28;
                rowIdx = 4;
                break;
            case 65:
                colIdx = 29;
                rowIdx = 4;
                break;
            case 66:
                colIdx = 30;
                rowIdx = 4;
                break;
            case 67:
                colIdx = 31;
                rowIdx = 4;
                break;

            default:
                colIdx = -99;
                rowIdx = -99;
                return;
        }
    } else if (agetIdx == 2) {
        switch (chanIdx) {
            case 0:
                colIdx = 8;
                rowIdx = 7;
                break;
            case 1:
                colIdx = 8;
                rowIdx = 6;
                break;
            case 2:
                colIdx = 8;
                rowIdx = 5;
                break;
            case 3:
                colIdx = 8;
                rowIdx = 4;
                break;
            case 4:
                colIdx = 8;
                rowIdx = 3;
                break;
            case 5:
                colIdx = 8;
                rowIdx = 2;
                break;
            case 6:
                colIdx = 8;
                rowIdx = 1;
                break;
            case 7:
                colIdx = 8;
                rowIdx = 0;
                break;
            case 8:
                colIdx = 9;
                rowIdx = 7;
                break;
            case 9:
                colIdx = 9;
                rowIdx = 6;
                break;
            case 10:
                colIdx = 9;
                rowIdx = 5;
                break;
            case 12:
                colIdx = 9;
                rowIdx = 4;
                break;
            case 13:
                colIdx = 9;
                rowIdx = 3;
                break;
            case 14:
                colIdx = 9;
                rowIdx = 2;
                break;
            case 15:
                colIdx = 9;
                rowIdx = 1;
                break;
            case 16:
                colIdx = 9;
                rowIdx = 0;
                break;
            case 17:
                colIdx = 10;
                rowIdx = 7;
                break;
            case 18:
                colIdx = 10;
                rowIdx = 6;
                break;
            case 19:
                colIdx = 10;
                rowIdx = 5;
                break;
            case 20:
                colIdx = 10;
                rowIdx = 4;
                break;
            case 21:
                colIdx = 10;
                rowIdx = 3;
                break;
            case 23:
                colIdx = 10;
                rowIdx = 2;
                break;
            case 24:
                colIdx = 10;
                rowIdx = 1;
                break;
            case 25:
                colIdx = 10;
                rowIdx = 0;
                break;
            case 26:
                colIdx = 11;
                rowIdx = 7;
                break;
            case 27:
                colIdx = 11;
                rowIdx = 6;
                break;
            case 28:
                colIdx = 11;
                rowIdx = 5;
                break;
            case 29:
                colIdx = 11;
                rowIdx = 4;
                break;
            case 30:
                colIdx = 11;
                rowIdx = 3;
                break;
            case 31:
                colIdx = 11;
                rowIdx = 2;
                break;
            case 32:
                colIdx = 11;
                rowIdx = 1;
                break;
            case 33:
                colIdx = 11;
                rowIdx = 0;
                break;
            case 34:
                colIdx = 12;
                rowIdx = 7;
                break;
            case 35:
                colIdx = 12;
                rowIdx = 6;
                break;
            case 36:
                colIdx = 12;
                rowIdx = 5;
                break;
            case 37:
                colIdx = 12;
                rowIdx = 4;
                break;
            case 38:
                colIdx = 12;
                rowIdx = 3;
                break;
            case 39:
                colIdx = 12;
                rowIdx = 2;
                break;
            case 40:
                colIdx = 12;
                rowIdx = 1;
                break;
            case 41:
                colIdx = 12;
                rowIdx = 0;
                break;
            case 42:
                colIdx = 13;
                rowIdx = 7;
                break;
            case 43:
                colIdx = 13;
                rowIdx = 6;
                break;
            case 44:
                colIdx = 13;
                rowIdx = 5;
                break;
            case 46:
                colIdx = 13;
                rowIdx = 4;
                break;
            case 47:
                colIdx = 13;
                rowIdx = 3;
                break;
            case 48:
                colIdx = 13;
                rowIdx = 2;
                break;
            case 49:
                colIdx = 13;
                rowIdx = 1;
                break;
            case 50:
                colIdx = 13;
                rowIdx = 0;
                break;
            case 51:
                colIdx = 14;
                rowIdx = 7;
                break;
            case 52:
                colIdx = 14;
                rowIdx = 6;
                break;
            case 53:
                colIdx = 14;
                rowIdx = 5;
                break;
            case 54:
                colIdx = 14;
                rowIdx = 4;
                break;
            case 55:
                colIdx = 14;
                rowIdx = 3;
                break;
            case 57:
                colIdx = 14;
                rowIdx = 2;
                break;
            case 58:
                colIdx = 14;
                rowIdx = 1;
                break;
            case 59:
                colIdx = 14;
                rowIdx = 0;
                break;
            case 60:
                colIdx = 15;
                rowIdx = 7;
                break;
            case 61:
                colIdx = 15;
                rowIdx = 6;
                break;
            case 62:
                colIdx = 15;
                rowIdx = 5;
                break;
            case 63:
                colIdx = 15;
                rowIdx = 4;
                break;
            case 64:
                colIdx = 15;
                rowIdx = 3;
                break;
            case 65:
                colIdx = 15;
                rowIdx = 2;
                break;
            case 66:
                colIdx = 15;
                rowIdx = 1;
                break;
            case 67:
                colIdx = 15;
                rowIdx = 0;
                break;

            default:
                colIdx = -99;
                rowIdx = -99;
                return;
        }
    } else if (agetIdx == 3) {
        switch (chanIdx) {
            case 0:
                colIdx = 0;
                rowIdx = 7;
                break;
            case 1:
                colIdx = 0;
                rowIdx = 6;
                break;
            case 2:
                colIdx = 0;
                rowIdx = 5;
                break;
            case 3:
                colIdx = 0;
                rowIdx = 4;
                break;
            case 4:
                colIdx = 0;
                rowIdx = 3;
                break;
            case 5:
                colIdx = 0;
                rowIdx = 2;
                break;
            case 6:
                colIdx = 0;
                rowIdx = 1;
                break;
            case 7:
                colIdx = 0;
                rowIdx = 0;
                break;
            case 8:
                colIdx = 1;
                rowIdx = 7;
                break;
            case 9:
                colIdx = 1;
                rowIdx = 6;
                break;
            case 10:
                colIdx = 1;
                rowIdx = 5;
                break;
            case 12:
                colIdx = 1;
                rowIdx = 4;
                break;
            case 13:
                colIdx = 1;
                rowIdx = 3;
                break;
            case 14:
                colIdx = 1;
                rowIdx = 2;
                break;
            case 15:
                colIdx = 1;
                rowIdx = 1;
                break;
            case 16:
                colIdx = 1;
                rowIdx = 0;
                break;
            case 17:
                colIdx = 2;
                rowIdx = 7;
                break;
            case 18:
                colIdx = 2;
                rowIdx = 6;
                break;
            case 19:
                colIdx = 2;
                rowIdx = 5;
                break;
            case 20:
                colIdx = 2;
                rowIdx = 4;
                break;
            case 21:
                colIdx = 2;
                rowIdx = 3;
                break;
            case 23:
                colIdx = 2;
                rowIdx = 2;
                break;
            case 24:
                colIdx = 2;
                rowIdx = 1;
                break;
            case 25:
                colIdx = 2;
                rowIdx = 0;
                break;
            case 26:
                colIdx = 3;
                rowIdx = 7;
                break;
            case 27:
                colIdx = 3;
                rowIdx = 6;
                break;
            case 28:
                colIdx = 3;
                rowIdx = 5;
                break;
            case 29:
                colIdx = 3;
                rowIdx = 4;
                break;
            case 30:
                colIdx = 3;
                rowIdx = 3;
                break;
            case 31:
                colIdx = 3;
                rowIdx = 2;
                break;
            case 32:
                colIdx = 3;
                rowIdx = 1;
                break;
            case 33:
                colIdx = 3;
                rowIdx = 0;
                break;
            case 34:
                colIdx = 4;
                rowIdx = 7;
                break;
            case 35:
                colIdx = 4;
                rowIdx = 6;
                break;
            case 36:
                colIdx = 4;
                rowIdx = 5;
                break;
            case 37:
                colIdx = 4;
                rowIdx = 4;
                break;
            case 38:
                colIdx = 4;
                rowIdx = 3;
                break;
            case 39:
                colIdx = 4;
                rowIdx = 2;
                break;
            case 40:
                colIdx = 4;
                rowIdx = 1;
                break;
            case 41:
                colIdx = 4;
                rowIdx = 0;
                break;
            case 42:
                colIdx = 5;
                rowIdx = 7;
                break;
            case 43:
                colIdx = 5;
                rowIdx = 6;
                break;
            case 44:
                colIdx = 5;
                rowIdx = 5;
                break;
            case 46:
                colIdx = 5;
                rowIdx = 4;
                break;
            case 47:
                colIdx = 5;
                rowIdx = 3;
                break;
            case 48:
                colIdx = 5;
                rowIdx = 2;
                break;
            case 49:
                colIdx = 5;
                rowIdx = 1;
                break;
            case 50:
                colIdx = 5;
                rowIdx = 0;
                break;
            case 51:
                colIdx = 6;
                rowIdx = 7;
                break;
            case 52:
                colIdx = 6;
                rowIdx = 6;
                break;
            case 53:
                colIdx = 6;
                rowIdx = 5;
                break;
            case 54:
                colIdx = 6;
                rowIdx = 4;
                break;
            case 55:
                colIdx = 6;
                rowIdx = 3;
                break;
            case 57:
                colIdx = 6;
                rowIdx = 2;
                break;
            case 58:
                colIdx = 6;
                rowIdx = 1;
                break;
            case 59:
                colIdx = 6;
                rowIdx = 0;
                break;
            case 60:
                colIdx = 7;
                rowIdx = 7;
                break;
            case 61:
                colIdx = 7;
                rowIdx = 6;
                break;
            case 62:
                colIdx = 7;
                rowIdx = 5;
                break;
            case 63:
                colIdx = 7;
                rowIdx = 4;
                break;
            case 64:
                colIdx = 7;
                rowIdx = 3;
                break;
            case 65:
                colIdx = 7;
                rowIdx = 2;
                break;
            case 66:
                colIdx = 7;
                rowIdx = 1;
                break;
            case 67:
                colIdx = 7;
                rowIdx = 0;
                break;

            default:
                colIdx = -99;
                rowIdx = -99;
                return;
        }
    }
}
