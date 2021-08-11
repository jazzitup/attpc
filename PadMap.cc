#include "PadMap.h"

PadMap::PadMap() {
    int xId = 0;
    int yId = 0;
    for (int agetId = 0; agetId < 4; agetId++) {
        switch (agetId) {
            case 0:
                xId = 16;
                yId = 3;
                break;
            case 1:
                xId = 16;
                yId = 7;
                break;
            case 2:
                xId = 8;
                yId = 7;
                break;
            case 3:
                xId = 0;
                yId = 7;
                break;
        }
        for (int chanId = 0; chanId < 68; chanId++) {
            if (chanId == 11 || chanId == 22 || chanId == 45 || chanId == 56) {
                this->xId_[agetId][chanId] = -1;
                this->yId_[agetId][chanId] = -1;
                continue;
            }
            if (agetId == 0 || agetId == 1) {
                this->agetId_[xId][yId] = agetId;
                this->chanId_[xId][yId] = chanId;
                this->xId_[agetId][chanId] = xId++;
                this->yId_[agetId][chanId] = yId;
                if (xId == 32) {
                    xId = 16;
                    yId--;
                }
            } else if (agetId == 2 || agetId == 3) {
                this->agetId_[xId][yId] = agetId;
                this->chanId_[xId][yId] = chanId;
                this->xId_[agetId][chanId] = xId;
                this->yId_[agetId][chanId] = yId--;
                if (yId == -1) {
                    xId++;
                    yId = 7;
                }
            }
        }
    }
}
int PadMap::GetAgetIdx(int xId, int yId) {
    return this->agetId_[xId][yId];
}
int PadMap::GetChanIdx(int xId, int yId) {
    return this->chanId_[xId][yId];
}
int PadMap::GetXId(int agetId, int chanId) {
    return this->xId_[agetId][chanId];
}
int PadMap::GetYId(int agetId, int chanId) {
    return this->yId_[agetId][chanId];
}
double PadMap::GetX(int agetId, int chanId) {
    return (this->dx_ + this->gap_) * this->GetXId(agetId, chanId);
}
double PadMap::GetY(int agetId, int chanId) {
    return (this->dy_ + this->gap_) * this->GetYId(agetId, chanId);
}
void PadMap::BuildPad(TH2Poly *poly) {
    // Add the bins
    double x[4], y[4];
    // Set Rectangle pattern parameters
    double xCenter = 0., yCenter = 0.;
    for (int xId = 0; xId < 32; xId++) {
        for (int yId = 0; yId < 8; yId++) {
            x[0] = xCenter - 0.5 * this->dx_;
            y[0] = yCenter - 0.5 * this->dy_;
            x[1] = x[0];
            y[1] = yCenter + 0.5 * this->dy_;
            x[2] = xCenter + 0.5 * this->dx_;
            y[2] = y[1];
            x[3] = x[2];
            y[3] = y[0];
            poly->AddBin(4, x, y);
            // Go up
            yCenter += this->dy_ + this->gap_;
        }
        xCenter += this->dx_ + this->gap_;
        yCenter = 0.;
    }
}
bool PadMap::IsDeadPad(int xId, int yId) {
    return this->IsDeadChan(agetId_[xId][yId], chanId_[xId][yId]);
}
bool PadMap::IsEvenChan(int agetId, int chanId) {
    if (chanId == 11 || chanId == 22 || chanId == 45 || chanId == 56 || chanId < 0 || chanId > 67 || agetId < 0 || agetId > 3) {
        std::cerr << "! Error in GETAnalyzer::IsEvenChan" << std::endl;
        std::cerr << "Out of range of agetId or chanId!" << std::endl;
        exit(-1);
    }
    if (chanId < 11 || (22 < chanId && chanId < 45) || 56 < chanId) {
        return (agetId == 0) ? (chanId % 2 == 0) : (chanId % 2 == 1);
    } else {
        return (agetId == 0) ? (chanId % 2 == 1) : (chanId % 2 == 0);
    }
}
bool PadMap::IsDeadChan(int agetId, int chanId) {
    if (chanId < 0 || chanId > 67 || agetId < 0 || agetId > 3) {
        std::cerr << "! Error in PadMap::IsDeadChan" << std::endl;
        std::cerr << "Out of range of agetId or chanId!" << std::endl;
        exit(-1);
    }
    if (chanId == 11 || chanId == 22 || chanId == 45 || chanId == 56) return false;
    // 52 (15 evens + 37 odds) channels is dead.
    if (agetId == 0) {
        switch (chanId) {  // AGET 0 even channels : 5 channels is dead.
            case 17:
                return true;
            case 32:
                return true;
            case 55:
                return true;
            case 58:
                return true;
            case 60:
                return true;
            default:
                return false;
        }
    } else if (agetId == 1) {
        if (this->IsEvenChan(agetId, chanId)) {
            switch (chanId) {  // AGET 1 even channels : 9 channels is dead.
                case 3:
                    return true;
                case 5:
                    return true;
                case 7:
                    return true;
                case 9:
                    return true;
                case 12:
                    return true;
                case 14:
                    return true;
                case 16:
                    return true;
                case 18:
                    return true;
                case 23:
                    return true;
                default:
                    return false;
            }
        } else {
            switch (chanId) {  // AGET 1 odd channels : 22 channels is dead.
                case 0:
                    return true;
                case 2:
                    return true;
                case 8:
                    return true;
                case 13:
                    return true;
                case 15:
                    return true;
                case 17:
                    return true;
                case 21:
                    return true;
                case 24:
                    return true;
                case 26:
                    return true;
                case 28:
                    return true;
                case 30:
                    return true;
                case 32:
                    return true;
                case 34:
                    return true;
                case 36:
                    return true;
                case 38:
                    return true;
                case 40:
                    return true;
                case 44:
                    return true;
                case 47:
                    return true;
                case 49:
                    return true;
                case 51:
                    return true;
                case 55:
                    return true;
                case 58:
                    return true;
                default:
                    return false;
            }
        }
    } else if (agetId == 2) {  // AGET 2 : all channels are alive.
        return false;
    } else {
        if (this->IsEvenChan(agetId, chanId)) {  // AGET 3 even channels : 3 channels is dead.
            switch (chanId) {
                case 39:
                    return true;
                case 61:
                    return true;
                case 65:
                    return true;
                default:
                    return false;
            }
        } else {
            switch (chanId) {  // AGET 3 odd channels : 15 channels is dead.
                case 0:
                    return true;
                case 2:
                    return true;
                case 6:
                    return true;
                case 8:
                    return true;
                case 10:
                    return true;
                case 13:
                    return true;
                case 15:
                    return true;
                case 17:
                    return true;
                case 21:
                    return true;
                case 24:
                    return true;
                case 26:
                    return true;
                case 28:
                    return true;
                case 34:
                    return true;
                case 38:
                    return true;
                case 40:
                    return true;
                case 42:
                    return true;
                default:
                    return false;
            }
        }
    }
}