#include "PadMap.h"
//
PadMap::PadMap(int padID) {
    this->padID_ = padID;
    std::vector<int> IDs;
    std::vector<int> position;
    int agetID, channelID;
    int colN, rowN;
    for (agetID = 0; agetID < 4; agetID++) {
        for (channelID = 0; channelID < 68; channelID++) {
            if (channelID == 11 || channelID == 22 || channelID == 45 || channelID == 56) continue;  // FPN Channel
            TransformChannelToPad(agetID, channelID, colN, rowN);
            IDs.clear();
            IDs.push_back(agetID);
            IDs.push_back(channelID);
            position.clear();
            position.push_back(colN);
            position.push_back(rowN);
            this->mapToID_[std::make_pair(colN, rowN)] = IDs;
            this->mapToPosition_[std::make_pair(agetID, channelID)] = position;
        }
    }
}
int PadMap::GetAgetIdx(int colN, int rowN) {
    return mapToID_[std::make_pair(colN, rowN)].at(0);
}
int PadMap::GetChanIdx(int colN, int rowN) {
    return mapToID_[std::make_pair(colN, rowN)].at(1);
}
int PadMap::GetColN(int agetID, int channelID) {
    return mapToPosition_[std::make_pair(agetID, channelID)].at(0);
}
int PadMap::GetRowN(int agetID, int channelID) {
    return mapToPosition_[std::make_pair(agetID, channelID)].at(1);
}
float PadMap::GetX(int colN) {
    if (this->padID_ == 0) {
        return 6.2 * colN;
    } else if (this->padID_ == 1) {
        return 3.125 * colN;
    }
    return -1;
}
float PadMap::GetY(int colN, int rowN) {
    if (this->padID_ == 0) {
        return colN % 2 ? 6.0 * rowN - 3.0 : 6.0 * rowN;
    } else if (this->padID_ == 1) {
        return 12.5 * rowN;
    }
    return -1;
}
void PadMap::BuildPad(TH2Poly *poly) {
    if (this->padID_ == 0) {
        // Add the bins
        double x[6], y[6];
        // Set hexagon pattern parameters
        double par[3], delX[2];
        par[0] = 6.2;
        par[1] = 6.;
        par[2] = 0.5;  // x-spacing, y-spacing, gap spacing
        delX[0] = 0.5 * (1 - par[2] / par[1]) * (pow(par[0], 2) + pow(0.5 * par[1], 2)) / par[0];
        delX[1] = sqrt(pow(delX[0], 2) - pow(0.5 * par[1] - par[2], 2));
        double xCenter = 0, yCenter = 0;
        for (int nRow = 0; nRow < 16; nRow++) {
            //  xtemp = xloop; // Resets the temp variable
            for (int nColumn = 0; nColumn < 16; nColumn++) {
                // Go around the hexagon
                x[0] = xCenter - delX[1];
                y[0] = yCenter - 0.5 * (par[1] - par[2]);
                x[1] = xCenter - delX[0];
                y[1] = yCenter;
                x[2] = x[0];
                y[2] = yCenter + 0.5 * (par[1] - par[2]);
                x[3] = xCenter + delX[1];
                y[3] = y[2];
                x[4] = xCenter + delX[0];
                y[4] = y[1];
                x[5] = x[3];
                y[5] = y[0];

                poly->AddBin(6, x, y);

                // Go up
                yCenter += par[1];
            }
            // Increment the starting position
            xCenter += par[0];
            if (nRow % 2 == 0)
                yCenter = -0.5 * par[1];
            else
                yCenter = 0;
        }
    } else if (this->padID_ == 1) {
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
}
void PadMap::TransformChannelToPad(int &agetID, int &channelID, int &colN, int &rowN) {
    // HoneyComb Pad: 0
    // Rectangle Pad: 1
    if (this->padID_ == 0) {  // HoneyComb Pad
        if (agetID == 0) {
            switch (channelID) {
                case 0:
                    colN = 8;
                    rowN = 7;
                    return;
                case 1:
                    colN = 9;
                    rowN = 7;
                    return;
                case 2:
                    colN = 10;
                    rowN = 7;
                    return;
                case 3:
                    colN = 11;
                    rowN = 7;
                    return;
                case 4:
                    colN = 12;
                    rowN = 7;
                    return;
                case 5:
                    colN = 13;
                    rowN = 7;
                    return;
                case 6:
                    colN = 14;
                    rowN = 7;
                    return;
                case 7:
                    colN = 15;
                    rowN = 7;
                    return;
                case 8:
                    colN = 8;
                    rowN = 6;
                    return;
                case 9:
                    colN = 9;
                    rowN = 6;
                    return;
                case 10:
                    colN = 10;
                    rowN = 6;
                    return;
                case 12:
                    colN = 11;
                    rowN = 6;
                    return;
                case 13:
                    colN = 12;
                    rowN = 6;
                    return;
                case 14:
                    colN = 13;
                    rowN = 6;
                    return;
                case 15:
                    colN = 14;
                    rowN = 6;
                    return;
                case 16:
                    colN = 15;
                    rowN = 6;
                    return;
                case 17:
                    colN = 8;
                    rowN = 5;
                    return;
                case 18:
                    colN = 9;
                    rowN = 5;
                    return;
                case 19:
                    colN = 10;
                    rowN = 5;
                    return;
                case 20:
                    colN = 11;
                    rowN = 5;
                    return;
                case 21:
                    colN = 12;
                    rowN = 5;
                    return;
                case 23:
                    colN = 13;
                    rowN = 5;
                    return;
                case 24:
                    colN = 14;
                    rowN = 5;
                    return;
                case 25:
                    colN = 15;
                    rowN = 5;
                    return;
                case 26:
                    colN = 8;
                    rowN = 4;
                    return;
                case 27:
                    colN = 9;
                    rowN = 4;
                    return;
                case 28:
                    colN = 10;
                    rowN = 4;
                    return;
                case 29:
                    colN = 11;
                    rowN = 4;
                    return;
                case 30:
                    colN = 12;
                    rowN = 4;
                    return;
                case 31:
                    colN = 13;
                    rowN = 4;
                    return;
                case 32:
                    colN = 14;
                    rowN = 4;
                    return;
                case 33:
                    colN = 15;
                    rowN = 4;
                    return;
                case 34:
                    colN = 8;
                    rowN = 3;
                    return;
                case 35:
                    colN = 9;
                    rowN = 3;
                    return;
                case 36:
                    colN = 10;
                    rowN = 3;
                    return;
                case 37:
                    colN = 11;
                    rowN = 3;
                    return;
                case 38:
                    colN = 12;
                    rowN = 3;
                    return;
                case 39:
                    colN = 13;
                    rowN = 3;
                    return;
                case 40:
                    colN = 14;
                    rowN = 3;
                    return;
                case 41:
                    colN = 15;
                    rowN = 3;
                    return;
                case 42:
                    colN = 8;
                    rowN = 2;
                    return;
                case 43:
                    colN = 9;
                    rowN = 2;
                    return;
                case 44:
                    colN = 10;
                    rowN = 2;
                    return;
                case 46:
                    colN = 11;
                    rowN = 2;
                    return;
                case 47:
                    colN = 12;
                    rowN = 2;
                    return;
                case 48:
                    colN = 13;
                    rowN = 2;
                    return;
                case 49:
                    colN = 14;
                    rowN = 2;
                    return;
                case 50:
                    colN = 15;
                    rowN = 2;
                    return;
                case 51:
                    colN = 8;
                    rowN = 1;
                    return;
                case 52:
                    colN = 9;
                    rowN = 1;
                    return;
                case 53:
                    colN = 10;
                    rowN = 1;
                    return;
                case 54:
                    colN = 11;
                    rowN = 1;
                    return;
                case 55:
                    colN = 12;
                    rowN = 1;
                    return;
                case 57:
                    colN = 13;
                    rowN = 1;
                    return;
                case 58:
                    colN = 14;
                    rowN = 1;
                    return;
                case 59:
                    colN = 15;
                    rowN = 1;
                    return;
                case 60:
                    colN = 8;
                    rowN = 0;
                    return;
                case 61:
                    colN = 9;
                    rowN = 0;
                    return;
                case 62:
                    colN = 10;
                    rowN = 0;
                    return;
                case 63:
                    colN = 11;
                    rowN = 0;
                    return;
                case 64:
                    colN = 12;
                    rowN = 0;
                    return;
                case 65:
                    colN = 13;
                    rowN = 0;
                    return;
                case 66:
                    colN = 14;
                    rowN = 0;
                    return;
                case 67:
                    colN = 15;
                    rowN = 0;
                    return;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        } else if (agetID == 1) {
            switch (channelID) {
                case 0:
                    colN = 8;
                    rowN = 15;
                    return;
                case 1:
                    colN = 9;
                    rowN = 15;
                    return;
                case 2:
                    colN = 10;
                    rowN = 15;
                    return;
                case 3:
                    colN = 11;
                    rowN = 15;
                    return;
                case 4:
                    colN = 12;
                    rowN = 15;
                    return;
                case 5:
                    colN = 13;
                    rowN = 15;
                    return;
                case 6:
                    colN = 14;
                    rowN = 15;
                    return;
                case 7:
                    colN = 15;
                    rowN = 15;
                    return;
                case 8:
                    colN = 8;
                    rowN = 14;
                    return;
                case 9:
                    colN = 9;
                    rowN = 14;
                    return;
                case 10:
                    colN = 10;
                    rowN = 14;
                    return;
                case 12:
                    colN = 11;
                    rowN = 14;
                    return;
                case 13:
                    colN = 12;
                    rowN = 14;
                    return;
                case 14:
                    colN = 13;
                    rowN = 14;
                    return;
                case 15:
                    colN = 14;
                    rowN = 14;
                    return;
                case 16:
                    colN = 15;
                    rowN = 14;
                    return;
                case 17:
                    colN = 8;
                    rowN = 13;
                    return;
                case 18:
                    colN = 9;
                    rowN = 13;
                    return;
                case 19:
                    colN = 10;
                    rowN = 13;
                    return;
                case 20:
                    colN = 11;
                    rowN = 13;
                    return;
                case 21:
                    colN = 12;
                    rowN = 13;
                    return;
                case 23:
                    colN = 13;
                    rowN = 13;
                    return;
                case 24:
                    colN = 14;
                    rowN = 13;
                    return;
                case 25:
                    colN = 15;
                    rowN = 13;
                    return;
                case 26:
                    colN = 8;
                    rowN = 12;
                    return;
                case 27:
                    colN = 9;
                    rowN = 12;
                    return;
                case 28:
                    colN = 10;
                    rowN = 12;
                    return;
                case 29:
                    colN = 11;
                    rowN = 12;
                    return;
                case 30:
                    colN = 12;
                    rowN = 12;
                    return;
                case 31:
                    colN = 13;
                    rowN = 12;
                    return;
                case 32:
                    colN = 14;
                    rowN = 12;
                    return;
                case 33:
                    colN = 15;
                    rowN = 12;
                    return;
                case 34:
                    colN = 8;
                    rowN = 11;
                    return;
                case 35:
                    colN = 9;
                    rowN = 11;
                    return;
                case 36:
                    colN = 10;
                    rowN = 11;
                    return;
                case 37:
                    colN = 11;
                    rowN = 11;
                    return;
                case 38:
                    colN = 12;
                    rowN = 11;
                    return;
                case 39:
                    colN = 13;
                    rowN = 11;
                    return;
                case 40:
                    colN = 14;
                    rowN = 11;
                    return;
                case 41:
                    colN = 15;
                    rowN = 11;
                    return;
                case 42:
                    colN = 8;
                    rowN = 10;
                    return;
                case 43:
                    colN = 9;
                    rowN = 10;
                    return;
                case 44:
                    colN = 10;
                    rowN = 10;
                    return;
                case 46:
                    colN = 11;
                    rowN = 10;
                    return;
                case 47:
                    colN = 12;
                    rowN = 10;
                    return;
                case 48:
                    colN = 13;
                    rowN = 10;
                    return;
                case 49:
                    colN = 14;
                    rowN = 10;
                    return;
                case 50:
                    colN = 15;
                    rowN = 10;
                    return;
                case 51:
                    colN = 8;
                    rowN = 9;
                    return;
                case 52:
                    colN = 9;
                    rowN = 9;
                    return;
                case 53:
                    colN = 10;
                    rowN = 9;
                    return;
                case 54:
                    colN = 11;
                    rowN = 9;
                    return;
                case 55:
                    colN = 12;
                    rowN = 9;
                    return;
                case 57:
                    colN = 13;
                    rowN = 9;
                    return;
                case 58:
                    colN = 14;
                    rowN = 9;
                    return;
                case 59:
                    colN = 15;
                    rowN = 9;
                    return;
                case 60:
                    colN = 8;
                    rowN = 8;
                    return;
                case 61:
                    colN = 9;
                    rowN = 8;
                    return;
                case 62:
                    colN = 10;
                    rowN = 8;
                    return;
                case 63:
                    colN = 11;
                    rowN = 8;
                    return;
                case 64:
                    colN = 12;
                    rowN = 8;
                    return;
                case 65:
                    colN = 13;
                    rowN = 8;
                    return;
                case 66:
                    colN = 14;
                    rowN = 8;
                    return;
                case 67:
                    colN = 15;
                    rowN = 8;
                    return;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        } else if (agetID == 2) {
            switch (channelID) {
                case 0:
                    colN = 4;
                    rowN = 15;
                    return;
                case 1:
                    colN = 4;
                    rowN = 14;
                    return;
                case 2:
                    colN = 4;
                    rowN = 13;
                    return;
                case 3:
                    colN = 4;
                    rowN = 12;
                    return;
                case 4:
                    colN = 4;
                    rowN = 11;
                    return;
                case 5:
                    colN = 4;
                    rowN = 10;
                    return;
                case 6:
                    colN = 4;
                    rowN = 9;
                    return;
                case 7:
                    colN = 4;
                    rowN = 8;
                    return;
                case 8:
                    colN = 4;
                    rowN = 7;
                    return;
                case 9:
                    colN = 4;
                    rowN = 6;
                    return;
                case 10:
                    colN = 4;
                    rowN = 5;
                    return;
                case 12:
                    colN = 4;
                    rowN = 4;
                    return;
                case 13:
                    colN = 4;
                    rowN = 3;
                    return;
                case 14:
                    colN = 4;
                    rowN = 2;
                    return;
                case 15:
                    colN = 4;
                    rowN = 1;
                    return;
                case 16:
                    colN = 4;
                    rowN = 0;
                    return;
                case 17:
                    colN = 5;
                    rowN = 15;
                    return;
                case 18:
                    colN = 5;
                    rowN = 14;
                    return;
                case 19:
                    colN = 5;
                    rowN = 13;
                    return;
                case 20:
                    colN = 5;
                    rowN = 12;
                    return;
                case 21:
                    colN = 5;
                    rowN = 11;
                    return;
                case 23:
                    colN = 5;
                    rowN = 10;
                    return;
                case 24:
                    colN = 5;
                    rowN = 9;
                    return;
                case 25:
                    colN = 5;
                    rowN = 8;
                    return;
                case 26:
                    colN = 5;
                    rowN = 7;
                    return;
                case 27:
                    colN = 5;
                    rowN = 6;
                    return;
                case 28:
                    colN = 5;
                    rowN = 5;
                    return;
                case 29:
                    colN = 5;
                    rowN = 4;
                    return;
                case 30:
                    colN = 5;
                    rowN = 3;
                    return;
                case 31:
                    colN = 5;
                    rowN = 2;
                    return;
                case 32:
                    colN = 5;
                    rowN = 1;
                    return;
                case 33:
                    colN = 5;
                    rowN = 0;
                    return;
                case 34:
                    colN = 6;
                    rowN = 15;
                    return;
                case 35:
                    colN = 6;
                    rowN = 14;
                    return;
                case 36:
                    colN = 6;
                    rowN = 13;
                    return;
                case 37:
                    colN = 6;
                    rowN = 12;
                    return;
                case 38:
                    colN = 6;
                    rowN = 11;
                    return;
                case 39:
                    colN = 6;
                    rowN = 10;
                    return;
                case 40:
                    colN = 6;
                    rowN = 9;
                    return;
                case 41:
                    colN = 6;
                    rowN = 8;
                    return;
                case 42:
                    colN = 6;
                    rowN = 7;
                    return;
                case 43:
                    colN = 6;
                    rowN = 6;
                    return;
                case 44:
                    colN = 6;
                    rowN = 5;
                    return;
                case 46:
                    colN = 6;
                    rowN = 4;
                    return;
                case 47:
                    colN = 6;
                    rowN = 3;
                    return;
                case 48:
                    colN = 6;
                    rowN = 2;
                    return;
                case 49:
                    colN = 6;
                    rowN = 1;
                    return;
                case 50:
                    colN = 6;
                    rowN = 0;
                    return;
                case 51:
                    colN = 7;
                    rowN = 15;
                    return;
                case 52:
                    colN = 7;
                    rowN = 14;
                    return;
                case 53:
                    colN = 7;
                    rowN = 13;
                    return;
                case 54:
                    colN = 7;
                    rowN = 12;
                    return;
                case 55:
                    colN = 7;
                    rowN = 11;
                    return;
                case 57:
                    colN = 7;
                    rowN = 10;
                    return;
                case 58:
                    colN = 7;
                    rowN = 9;
                    return;
                case 59:
                    colN = 7;
                    rowN = 8;
                    return;
                case 60:
                    colN = 7;
                    rowN = 7;
                    return;
                case 61:
                    colN = 7;
                    rowN = 6;
                    return;
                case 62:
                    colN = 7;
                    rowN = 5;
                    return;
                case 63:
                    colN = 7;
                    rowN = 4;
                    return;
                case 64:
                    colN = 7;
                    rowN = 3;
                    return;
                case 65:
                    colN = 7;
                    rowN = 2;
                    return;
                case 66:
                    colN = 7;
                    rowN = 1;
                    return;
                case 67:
                    colN = 7;
                    rowN = 0;
                    return;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        } else if (agetID == 3) {
            switch (channelID) {
                case 0:
                    colN = 0;
                    rowN = 15;
                    return;
                case 1:
                    colN = 0;
                    rowN = 14;
                    return;
                case 2:
                    colN = 0;
                    rowN = 13;
                    return;
                case 3:
                    colN = 0;
                    rowN = 12;
                    return;
                case 4:
                    colN = 0;
                    rowN = 11;
                    return;
                case 5:
                    colN = 0;
                    rowN = 10;
                    return;
                case 6:
                    colN = 0;
                    rowN = 9;
                    return;
                case 7:
                    colN = 0;
                    rowN = 8;
                    return;
                case 8:
                    colN = 0;
                    rowN = 7;
                    return;
                case 9:
                    colN = 0;
                    rowN = 6;
                    return;
                case 10:
                    colN = 0;
                    rowN = 5;
                    return;
                case 12:
                    colN = 0;
                    rowN = 4;
                    return;
                case 13:
                    colN = 0;
                    rowN = 3;
                    return;
                case 14:
                    colN = 0;
                    rowN = 2;
                    return;
                case 15:
                    colN = 0;
                    rowN = 1;
                    return;
                case 16:
                    colN = 0;
                    rowN = 0;
                    return;
                case 17:
                    colN = 1;
                    rowN = 15;
                    return;
                case 18:
                    colN = 1;
                    rowN = 14;
                    return;
                case 19:
                    colN = 1;
                    rowN = 13;
                    return;
                case 20:
                    colN = 1;
                    rowN = 12;
                    return;
                case 21:
                    colN = 1;
                    rowN = 11;
                    return;
                case 23:
                    colN = 1;
                    rowN = 10;
                    return;
                case 24:
                    colN = 1;
                    rowN = 9;
                    return;
                case 25:
                    colN = 1;
                    rowN = 8;
                    return;
                case 26:
                    colN = 1;
                    rowN = 7;
                    return;
                case 27:
                    colN = 1;
                    rowN = 6;
                    return;
                case 28:
                    colN = 1;
                    rowN = 5;
                    return;
                case 29:
                    colN = 1;
                    rowN = 4;
                    return;
                case 30:
                    colN = 1;
                    rowN = 3;
                    return;
                case 31:
                    colN = 1;
                    rowN = 2;
                    return;
                case 32:
                    colN = 1;
                    rowN = 1;
                    return;
                case 33:
                    colN = 1;
                    rowN = 0;
                    return;
                case 34:
                    colN = 2;
                    rowN = 15;
                    return;
                case 35:
                    colN = 2;
                    rowN = 14;
                    return;
                case 36:
                    colN = 2;
                    rowN = 13;
                    return;
                case 37:
                    colN = 2;
                    rowN = 12;
                    return;
                case 38:
                    colN = 2;
                    rowN = 11;
                    return;
                case 39:
                    colN = 2;
                    rowN = 10;
                    return;
                case 40:
                    colN = 2;
                    rowN = 9;
                    return;
                case 41:
                    colN = 2;
                    rowN = 8;
                    return;
                case 42:
                    colN = 2;
                    rowN = 7;
                    return;
                case 43:
                    colN = 2;
                    rowN = 6;
                    return;
                case 44:
                    colN = 2;
                    rowN = 5;
                    return;
                case 46:
                    colN = 2;
                    rowN = 4;
                    return;
                case 47:
                    colN = 2;
                    rowN = 3;
                    return;
                case 48:
                    colN = 2;
                    rowN = 2;
                    return;
                case 49:
                    colN = 2;
                    rowN = 1;
                    return;
                case 50:
                    colN = 2;
                    rowN = 0;
                    return;
                case 51:
                    colN = 3;
                    rowN = 15;
                    return;
                case 52:
                    colN = 3;
                    rowN = 14;
                    return;
                case 53:
                    colN = 3;
                    rowN = 13;
                    return;
                case 54:
                    colN = 3;
                    rowN = 12;
                    return;
                case 55:
                    colN = 3;
                    rowN = 11;
                    return;
                case 57:
                    colN = 3;
                    rowN = 10;
                    return;
                case 58:
                    colN = 3;
                    rowN = 9;
                    return;
                case 59:
                    colN = 3;
                    rowN = 8;
                    return;
                case 60:
                    colN = 3;
                    rowN = 7;
                    return;
                case 61:
                    colN = 3;
                    rowN = 6;
                    return;
                case 62:
                    colN = 3;
                    rowN = 5;
                    return;
                case 63:
                    colN = 3;
                    rowN = 4;
                    return;
                case 64:
                    colN = 3;
                    rowN = 3;
                    return;
                case 65:
                    colN = 3;
                    rowN = 2;
                    return;
                case 66:
                    colN = 3;
                    rowN = 1;
                    return;
                case 67:
                    colN = 3;
                    rowN = 0;
                    return;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        }
    } else if (this->padID_ == 1) {  // Rectangle Pad
        if (agetID == 0) {
            switch (channelID) {
                case 0:
                    colN = 16;
                    rowN = 3;
                    break;
                case 1:
                    colN = 17;
                    rowN = 3;
                    break;
                case 2:
                    colN = 18;
                    rowN = 3;
                    break;
                case 3:
                    colN = 19;
                    rowN = 3;
                    break;
                case 4:
                    colN = 20;
                    rowN = 3;
                    break;
                case 5:
                    colN = 21;
                    rowN = 3;
                    break;
                case 6:
                    colN = 22;
                    rowN = 3;
                    break;
                case 7:
                    colN = 23;
                    rowN = 3;
                    break;
                case 8:
                    colN = 24;
                    rowN = 3;
                    break;
                case 9:
                    colN = 25;
                    rowN = 3;
                    break;
                case 10:
                    colN = 26;
                    rowN = 3;
                    break;
                case 12:
                    colN = 27;
                    rowN = 3;
                    break;
                case 13:
                    colN = 28;
                    rowN = 3;
                    break;
                case 14:
                    colN = 29;
                    rowN = 3;
                    break;
                case 15:
                    colN = 30;
                    rowN = 3;
                    break;
                case 16:
                    colN = 31;
                    rowN = 3;
                    break;
                case 17:
                    colN = 16;
                    rowN = 2;
                    break;
                case 18:
                    colN = 17;
                    rowN = 2;
                    break;
                case 19:
                    colN = 18;
                    rowN = 2;
                    break;
                case 20:
                    colN = 19;
                    rowN = 2;
                    break;
                case 21:
                    colN = 20;
                    rowN = 2;
                    break;
                case 23:
                    colN = 21;
                    rowN = 2;
                    break;
                case 24:
                    colN = 22;
                    rowN = 2;
                    break;
                case 25:
                    colN = 23;
                    rowN = 2;
                    break;
                case 26:
                    colN = 24;
                    rowN = 2;
                    break;
                case 27:
                    colN = 25;
                    rowN = 2;
                    break;
                case 28:
                    colN = 26;
                    rowN = 2;
                    break;
                case 29:
                    colN = 27;
                    rowN = 2;
                    break;
                case 30:
                    colN = 28;
                    rowN = 2;
                    break;
                case 31:
                    colN = 29;
                    rowN = 2;
                    break;
                case 32:
                    colN = 30;
                    rowN = 2;
                    break;
                case 33:
                    colN = 31;
                    rowN = 2;
                    break;
                case 34:
                    colN = 16;
                    rowN = 1;
                    break;
                case 35:
                    colN = 17;
                    rowN = 1;
                    break;
                case 36:
                    colN = 18;
                    rowN = 1;
                    break;
                case 37:
                    colN = 19;
                    rowN = 1;
                    break;
                case 38:
                    colN = 20;
                    rowN = 1;
                    break;
                case 39:
                    colN = 21;
                    rowN = 1;
                    break;
                case 40:
                    colN = 22;
                    rowN = 1;
                    break;
                case 41:
                    colN = 23;
                    rowN = 1;
                    break;
                case 42:
                    colN = 24;
                    rowN = 1;
                    break;
                case 43:
                    colN = 25;
                    rowN = 1;
                    break;
                case 44:
                    colN = 26;
                    rowN = 1;
                    break;
                case 46:
                    colN = 27;
                    rowN = 1;
                    break;
                case 47:
                    colN = 28;
                    rowN = 1;
                    break;
                case 48:
                    colN = 29;
                    rowN = 1;
                    break;
                case 49:
                    colN = 30;
                    rowN = 1;
                    break;
                case 50:
                    colN = 31;
                    rowN = 1;
                    break;
                case 51:
                    colN = 16;
                    rowN = 0;
                    break;
                case 52:
                    colN = 17;
                    rowN = 0;
                    break;
                case 53:
                    colN = 18;
                    rowN = 0;
                    break;
                case 54:
                    colN = 19;
                    rowN = 0;
                    break;
                case 55:
                    colN = 20;
                    rowN = 0;
                    break;
                case 57:
                    colN = 21;
                    rowN = 0;
                    break;
                case 58:
                    colN = 22;
                    rowN = 0;
                    break;
                case 59:
                    colN = 23;
                    rowN = 0;
                    break;
                case 60:
                    colN = 24;
                    rowN = 0;
                    break;
                case 61:
                    colN = 25;
                    rowN = 0;
                    break;
                case 62:
                    colN = 26;
                    rowN = 0;
                    break;
                case 63:
                    colN = 27;
                    rowN = 0;
                    break;
                case 64:
                    colN = 28;
                    rowN = 0;
                    break;
                case 65:
                    colN = 29;
                    rowN = 0;
                    break;
                case 66:
                    colN = 30;
                    rowN = 0;
                    break;
                case 67:
                    colN = 31;
                    rowN = 0;
                    break;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        } else if (agetID == 1) {
            switch (channelID) {
                case 0:
                    colN = 16;
                    rowN = 7;
                    break;
                case 1:
                    colN = 17;
                    rowN = 7;
                    break;
                case 2:
                    colN = 18;
                    rowN = 7;
                    break;
                case 3:
                    colN = 19;
                    rowN = 7;
                    break;
                case 4:
                    colN = 20;
                    rowN = 7;
                    break;
                case 5:
                    colN = 21;
                    rowN = 7;
                    break;
                case 6:
                    colN = 22;
                    rowN = 7;
                    break;
                case 7:
                    colN = 23;
                    rowN = 7;
                    break;
                case 8:
                    colN = 24;
                    rowN = 7;
                    break;
                case 9:
                    colN = 25;
                    rowN = 7;
                    break;
                case 10:
                    colN = 26;
                    rowN = 7;
                    break;
                case 12:
                    colN = 27;
                    rowN = 7;
                    break;
                case 13:
                    colN = 28;
                    rowN = 7;
                    break;
                case 14:
                    colN = 29;
                    rowN = 7;
                    break;
                case 15:
                    colN = 30;
                    rowN = 7;
                    break;
                case 16:
                    colN = 31;
                    rowN = 7;
                    break;
                case 17:
                    colN = 16;
                    rowN = 6;
                    break;
                case 18:
                    colN = 17;
                    rowN = 6;
                    break;
                case 19:
                    colN = 18;
                    rowN = 6;
                    break;
                case 20:
                    colN = 19;
                    rowN = 6;
                    break;
                case 21:
                    colN = 20;
                    rowN = 6;
                    break;
                case 23:
                    colN = 21;
                    rowN = 6;
                    break;
                case 24:
                    colN = 22;
                    rowN = 6;
                    break;
                case 25:
                    colN = 23;
                    rowN = 6;
                    break;
                case 26:
                    colN = 24;
                    rowN = 6;
                    break;
                case 27:
                    colN = 25;
                    rowN = 6;
                    break;
                case 28:
                    colN = 26;
                    rowN = 6;
                    break;
                case 29:
                    colN = 27;
                    rowN = 6;
                    break;
                case 30:
                    colN = 28;
                    rowN = 6;
                    break;
                case 31:
                    colN = 29;
                    rowN = 6;
                    break;
                case 32:
                    colN = 30;
                    rowN = 6;
                    break;
                case 33:
                    colN = 31;
                    rowN = 6;
                    break;
                case 34:
                    colN = 16;
                    rowN = 5;
                    break;
                case 35:
                    colN = 17;
                    rowN = 5;
                    break;
                case 36:
                    colN = 18;
                    rowN = 5;
                    break;
                case 37:
                    colN = 19;
                    rowN = 5;
                    break;
                case 38:
                    colN = 20;
                    rowN = 5;
                    break;
                case 39:
                    colN = 21;
                    rowN = 5;
                    break;
                case 40:
                    colN = 22;
                    rowN = 5;
                    break;
                case 41:
                    colN = 23;
                    rowN = 5;
                    break;
                case 42:
                    colN = 24;
                    rowN = 5;
                    break;
                case 43:
                    colN = 25;
                    rowN = 5;
                    break;
                case 44:
                    colN = 26;
                    rowN = 5;
                    break;
                case 46:
                    colN = 27;
                    rowN = 5;
                    break;
                case 47:
                    colN = 28;
                    rowN = 5;
                    break;
                case 48:
                    colN = 29;
                    rowN = 5;
                    break;
                case 49:
                    colN = 30;
                    rowN = 5;
                    break;
                case 50:
                    colN = 31;
                    rowN = 5;
                    break;
                case 51:
                    colN = 16;
                    rowN = 4;
                    break;
                case 52:
                    colN = 17;
                    rowN = 4;
                    break;
                case 53:
                    colN = 18;
                    rowN = 4;
                    break;
                case 54:
                    colN = 19;
                    rowN = 4;
                    break;
                case 55:
                    colN = 20;
                    rowN = 4;
                    break;
                case 57:
                    colN = 21;
                    rowN = 4;
                    break;
                case 58:
                    colN = 22;
                    rowN = 4;
                    break;
                case 59:
                    colN = 23;
                    rowN = 4;
                    break;
                case 60:
                    colN = 24;
                    rowN = 4;
                    break;
                case 61:
                    colN = 25;
                    rowN = 4;
                    break;
                case 62:
                    colN = 26;
                    rowN = 4;
                    break;
                case 63:
                    colN = 27;
                    rowN = 4;
                    break;
                case 64:
                    colN = 28;
                    rowN = 4;
                    break;
                case 65:
                    colN = 29;
                    rowN = 4;
                    break;
                case 66:
                    colN = 30;
                    rowN = 4;
                    break;
                case 67:
                    colN = 31;
                    rowN = 4;
                    break;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        } else if (agetID == 2) {
            switch (channelID) {
                case 0:
                    colN = 8;
                    rowN = 7;
                    break;
                case 1:
                    colN = 8;
                    rowN = 6;
                    break;
                case 2:
                    colN = 8;
                    rowN = 5;
                    break;
                case 3:
                    colN = 8;
                    rowN = 4;
                    break;
                case 4:
                    colN = 8;
                    rowN = 3;
                    break;
                case 5:
                    colN = 8;
                    rowN = 2;
                    break;
                case 6:
                    colN = 8;
                    rowN = 1;
                    break;
                case 7:
                    colN = 8;
                    rowN = 0;
                    break;
                case 8:
                    colN = 9;
                    rowN = 7;
                    break;
                case 9:
                    colN = 9;
                    rowN = 6;
                    break;
                case 10:
                    colN = 9;
                    rowN = 5;
                    break;
                case 12:
                    colN = 9;
                    rowN = 4;
                    break;
                case 13:
                    colN = 9;
                    rowN = 3;
                    break;
                case 14:
                    colN = 9;
                    rowN = 2;
                    break;
                case 15:
                    colN = 9;
                    rowN = 1;
                    break;
                case 16:
                    colN = 9;
                    rowN = 0;
                    break;
                case 17:
                    colN = 10;
                    rowN = 7;
                    break;
                case 18:
                    colN = 10;
                    rowN = 6;
                    break;
                case 19:
                    colN = 10;
                    rowN = 5;
                    break;
                case 20:
                    colN = 10;
                    rowN = 4;
                    break;
                case 21:
                    colN = 10;
                    rowN = 3;
                    break;
                case 23:
                    colN = 10;
                    rowN = 2;
                    break;
                case 24:
                    colN = 10;
                    rowN = 1;
                    break;
                case 25:
                    colN = 10;
                    rowN = 0;
                    break;
                case 26:
                    colN = 11;
                    rowN = 7;
                    break;
                case 27:
                    colN = 11;
                    rowN = 6;
                    break;
                case 28:
                    colN = 11;
                    rowN = 5;
                    break;
                case 29:
                    colN = 11;
                    rowN = 4;
                    break;
                case 30:
                    colN = 11;
                    rowN = 3;
                    break;
                case 31:
                    colN = 11;
                    rowN = 2;
                    break;
                case 32:
                    colN = 11;
                    rowN = 1;
                    break;
                case 33:
                    colN = 11;
                    rowN = 0;
                    break;
                case 34:
                    colN = 12;
                    rowN = 7;
                    break;
                case 35:
                    colN = 12;
                    rowN = 6;
                    break;
                case 36:
                    colN = 12;
                    rowN = 5;
                    break;
                case 37:
                    colN = 12;
                    rowN = 4;
                    break;
                case 38:
                    colN = 12;
                    rowN = 3;
                    break;
                case 39:
                    colN = 12;
                    rowN = 2;
                    break;
                case 40:
                    colN = 12;
                    rowN = 1;
                    break;
                case 41:
                    colN = 12;
                    rowN = 0;
                    break;
                case 42:
                    colN = 13;
                    rowN = 7;
                    break;
                case 43:
                    colN = 13;
                    rowN = 6;
                    break;
                case 44:
                    colN = 13;
                    rowN = 5;
                    break;
                case 46:
                    colN = 13;
                    rowN = 4;
                    break;
                case 47:
                    colN = 13;
                    rowN = 3;
                    break;
                case 48:
                    colN = 13;
                    rowN = 2;
                    break;
                case 49:
                    colN = 13;
                    rowN = 1;
                    break;
                case 50:
                    colN = 13;
                    rowN = 0;
                    break;
                case 51:
                    colN = 14;
                    rowN = 7;
                    break;
                case 52:
                    colN = 14;
                    rowN = 6;
                    break;
                case 53:
                    colN = 14;
                    rowN = 5;
                    break;
                case 54:
                    colN = 14;
                    rowN = 4;
                    break;
                case 55:
                    colN = 14;
                    rowN = 3;
                    break;
                case 57:
                    colN = 14;
                    rowN = 2;
                    break;
                case 58:
                    colN = 14;
                    rowN = 1;
                    break;
                case 59:
                    colN = 14;
                    rowN = 0;
                    break;
                case 60:
                    colN = 15;
                    rowN = 7;
                    break;
                case 61:
                    colN = 15;
                    rowN = 6;
                    break;
                case 62:
                    colN = 15;
                    rowN = 5;
                    break;
                case 63:
                    colN = 15;
                    rowN = 4;
                    break;
                case 64:
                    colN = 15;
                    rowN = 3;
                    break;
                case 65:
                    colN = 15;
                    rowN = 2;
                    break;
                case 66:
                    colN = 15;
                    rowN = 1;
                    break;
                case 67:
                    colN = 15;
                    rowN = 0;
                    break;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        } else if (agetID == 3) {
            switch (channelID) {
                case 0:
                    colN = 0;
                    rowN = 7;
                    break;
                case 1:
                    colN = 0;
                    rowN = 6;
                    break;
                case 2:
                    colN = 0;
                    rowN = 5;
                    break;
                case 3:
                    colN = 0;
                    rowN = 4;
                    break;
                case 4:
                    colN = 0;
                    rowN = 3;
                    break;
                case 5:
                    colN = 0;
                    rowN = 2;
                    break;
                case 6:
                    colN = 0;
                    rowN = 1;
                    break;
                case 7:
                    colN = 0;
                    rowN = 0;
                    break;
                case 8:
                    colN = 1;
                    rowN = 7;
                    break;
                case 9:
                    colN = 1;
                    rowN = 6;
                    break;
                case 10:
                    colN = 1;
                    rowN = 5;
                    break;
                case 12:
                    colN = 1;
                    rowN = 4;
                    break;
                case 13:
                    colN = 1;
                    rowN = 3;
                    break;
                case 14:
                    colN = 1;
                    rowN = 2;
                    break;
                case 15:
                    colN = 1;
                    rowN = 1;
                    break;
                case 16:
                    colN = 1;
                    rowN = 0;
                    break;
                case 17:
                    colN = 2;
                    rowN = 7;
                    break;
                case 18:
                    colN = 2;
                    rowN = 6;
                    break;
                case 19:
                    colN = 2;
                    rowN = 5;
                    break;
                case 20:
                    colN = 2;
                    rowN = 4;
                    break;
                case 21:
                    colN = 2;
                    rowN = 3;
                    break;
                case 23:
                    colN = 2;
                    rowN = 2;
                    break;
                case 24:
                    colN = 2;
                    rowN = 1;
                    break;
                case 25:
                    colN = 2;
                    rowN = 0;
                    break;
                case 26:
                    colN = 3;
                    rowN = 7;
                    break;
                case 27:
                    colN = 3;
                    rowN = 6;
                    break;
                case 28:
                    colN = 3;
                    rowN = 5;
                    break;
                case 29:
                    colN = 3;
                    rowN = 4;
                    break;
                case 30:
                    colN = 3;
                    rowN = 3;
                    break;
                case 31:
                    colN = 3;
                    rowN = 2;
                    break;
                case 32:
                    colN = 3;
                    rowN = 1;
                    break;
                case 33:
                    colN = 3;
                    rowN = 0;
                    break;
                case 34:
                    colN = 4;
                    rowN = 7;
                    break;
                case 35:
                    colN = 4;
                    rowN = 6;
                    break;
                case 36:
                    colN = 4;
                    rowN = 5;
                    break;
                case 37:
                    colN = 4;
                    rowN = 4;
                    break;
                case 38:
                    colN = 4;
                    rowN = 3;
                    break;
                case 39:
                    colN = 4;
                    rowN = 2;
                    break;
                case 40:
                    colN = 4;
                    rowN = 1;
                    break;
                case 41:
                    colN = 4;
                    rowN = 0;
                    break;
                case 42:
                    colN = 5;
                    rowN = 7;
                    break;
                case 43:
                    colN = 5;
                    rowN = 6;
                    break;
                case 44:
                    colN = 5;
                    rowN = 5;
                    break;
                case 46:
                    colN = 5;
                    rowN = 4;
                    break;
                case 47:
                    colN = 5;
                    rowN = 3;
                    break;
                case 48:
                    colN = 5;
                    rowN = 2;
                    break;
                case 49:
                    colN = 5;
                    rowN = 1;
                    break;
                case 50:
                    colN = 5;
                    rowN = 0;
                    break;
                case 51:
                    colN = 6;
                    rowN = 7;
                    break;
                case 52:
                    colN = 6;
                    rowN = 6;
                    break;
                case 53:
                    colN = 6;
                    rowN = 5;
                    break;
                case 54:
                    colN = 6;
                    rowN = 4;
                    break;
                case 55:
                    colN = 6;
                    rowN = 3;
                    break;
                case 57:
                    colN = 6;
                    rowN = 2;
                    break;
                case 58:
                    colN = 6;
                    rowN = 1;
                    break;
                case 59:
                    colN = 6;
                    rowN = 0;
                    break;
                case 60:
                    colN = 7;
                    rowN = 7;
                    break;
                case 61:
                    colN = 7;
                    rowN = 6;
                    break;
                case 62:
                    colN = 7;
                    rowN = 5;
                    break;
                case 63:
                    colN = 7;
                    rowN = 4;
                    break;
                case 64:
                    colN = 7;
                    rowN = 3;
                    break;
                case 65:
                    colN = 7;
                    rowN = 2;
                    break;
                case 66:
                    colN = 7;
                    rowN = 1;
                    break;
                case 67:
                    colN = 7;
                    rowN = 0;
                    break;

                default:
                    colN = -99;
                    rowN = -99;
                    return;
            }
        }
    }
}