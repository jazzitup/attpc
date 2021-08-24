#include "GETAnalyzer.h"

GETAnalyzer::GETAnalyzer() {
}
void GETAnalyzer::LinkToDecoder(GETDecoder *decoder) {
    this->decoder_ = decoder;
}
void GETAnalyzer::RunEventChecker() {
    if (this->decoder_->GetEventId() == 0) {
        this->nFakeEvent_ = 0;
        this->nSparkEvent_ = 0;
    }
    // *** Check whether it's fake or not. ***
    // if 10 ms < diffTime < 30 ms, this event is a fake trigger event.
    this->isFakeEvent_ = false;
    if (1E6 < this->decoder_->GetDiffTime() && this->decoder_->GetDiffTime() < 3E6) {
        this->isFakeEvent_ = true;
        this->nFakeEvent_++;
    }
    if (this->decoder_->GetEventId() == 0) this->isFakeEvent_ = false;
    // *** Check whether it's spark or not. ***
    // if At least one thing is mimimum ADC < 100, this event is a spark event.
    this->isSparkEvent_ = false;
    for (int agetId = 0; agetId < 4; agetId++) {
        for (int chanId = 0; chanId < 68; chanId++) {
            if (this->IsFPNChan(chanId) || this->IsDeadChan(agetId, chanId)) continue;
            if (this->decoder_->GetMinADC(agetId, chanId) < 100) {
                this->isSparkEvent_ = true;
                this->nSparkEvent_++;
                return;
            }
        }
    }
}
void GETAnalyzer::RunNoiseCanceller() {
    // Step 1: Subtract a Fixed Pattern Noise from raw signals
    for (int agetId = 0; agetId < 4; agetId++) {
        int nChanId = 0;
        for (int chanId = 0; chanId < 68; chanId++) {
            if (this->IsFPNChan(chanId)) continue;
            nChanId++;
            int FPN_id = (nChanId - 1) / 16;
            if (this->IsDeadChan(agetId, chanId)) continue;
            for (int buckId = 0; buckId < 512; buckId++) {
                double ADC_scaled = this->decoder_->GetADC(agetId, chanId, buckId);
                double FPN_scaled = this->decoder_->GetADC(agetId, FPN_chanId_[FPN_id], buckId);
                this->ADC_FPNC_[agetId][chanId][buckId] = ADC_scaled - FPN_scaled;
            }
        }
    }
    // Step 2: Subtract a mean value from the (FPNC)signals to find noise candidates
    for (int agetId = 0; agetId < 4; agetId++) {
        for (int chanId = 0; chanId < 68; chanId++) {
            if (this->IsFPNChan(chanId) || this->IsDeadChan(agetId, chanId)) continue;
            double mean = 0.;
            for (int buckId = 1; buckId <= 500; buckId++) mean += this->ADC_FPNC_[agetId][chanId][buckId] / 500.;
            for (int buckId = 0; buckId < 512; buckId++) this->ADC_FPNC_[agetId][chanId][buckId] -= mean;
        }
    }
    // Step 3: make noise profiles using 4 noise candidates
    for (int agetId = 0; agetId < 4; agetId++) {
        std::vector<std::pair<int, double>> vChanId_maxADC_evn;
        std::vector<std::pair<int, double>> vChanId_maxADC_odd;
        for (int chanId = 0; chanId < 68; chanId++) {
            if (this->IsFPNChan(chanId) || this->IsDeadChan(agetId, chanId)) continue;
            if (this->IsEvenChan(agetId, chanId))
                vChanId_maxADC_evn.push_back(std::make_pair(chanId, this->ADC_FPNC_[agetId][chanId][1]));
            else
                vChanId_maxADC_odd.push_back(std::make_pair(chanId, this->ADC_FPNC_[agetId][chanId][1]));
        }
        // Sort by a descending order
        sort(vChanId_maxADC_evn.begin(), vChanId_maxADC_evn.end(), compareSecondByDescending);
        sort(vChanId_maxADC_odd.begin(), vChanId_maxADC_odd.end(), compareSecondByDescending);
        for (int buckId = 0; buckId < 512; buckId++) {
            this->noiseEvn_[agetId][buckId] = 0.;
            this->noiseOdd_[agetId][buckId] = 0.;
            for (int idx = 0; idx < 4; idx++) {
                int evnChanId = vChanId_maxADC_evn[idx].first;
                int oddChanId = vChanId_maxADC_odd[idx].first;
                this->noiseEvn_[agetId][buckId] += 0.25 * this->ADC_FPNC_[agetId][evnChanId][buckId];
                this->noiseOdd_[agetId][buckId] += 0.25 * this->ADC_FPNC_[agetId][oddChanId][buckId];
            }
        }
    }
    // Step 4: adjust a pedestal of the (FPNC)signals to 0
    for (int agetId = 0; agetId < 4; agetId++) {
        for (int chanId = 0; chanId < 68; chanId++) {
            if (this->IsFPNChan(chanId) || this->IsDeadChan(agetId, chanId)) continue;
            std::vector<double> vADC;
            for (int buckId = 1; buckId <= 500; buckId++) vADC.push_back(this->ADC_FPNC_[agetId][chanId][buckId]);
            std::sort(vADC.begin(), vADC.end());
            double deviation = 0.;
            for (int buckId = 151; buckId <= 200; buckId++) deviation += vADC[buckId - 1] / 50.;
            for (int buckId = 0; buckId < 512; buckId++) this->ADC_FPNC_[agetId][chanId][buckId] -= deviation;
        }
    }
    // Step 5: subtract a noise profile from the (FPNC)signals
    for (int agetId = 0; agetId < 4; agetId++) {
        for (int chanId = 0; chanId < 68; chanId++) {
            if (this->IsFPNChan(chanId) || this->IsDeadChan(agetId, chanId)) continue;
            for (int buckId = 0; buckId < 512; buckId++) {
                this->ADC_ANC_[agetId][chanId][buckId] = this->ADC_FPNC_[agetId][chanId][buckId];
                if (this->IsEvenChan(agetId, chanId)) {
                    this->ADC_ANC_[agetId][chanId][buckId] -= this->noiseEvn_[agetId][buckId];
                } else {
                    this->ADC_ANC_[agetId][chanId][buckId] -= this->noiseOdd_[agetId][buckId];
                }
            }
        }
    }
}
double GETAnalyzer::GetADC_ANC(int agetId, int chanId, int buckId) {
    return this->ADC_ANC_[agetId][chanId][buckId];
}
double GETAnalyzer::GetADC_FPNC(int agetId, int chanId, int buckId) {
    return this->ADC_FPNC_[agetId][chanId][buckId];
}
double GETAnalyzer::GetNoise(int agetId, int buckId, bool useEvn) {
    return useEvn ? this->noiseEvn_[agetId][buckId] : this->noiseOdd_[agetId][buckId];
}
bool GETAnalyzer::IsFakeEvent() {
    return this->isFakeEvent_;
}
bool GETAnalyzer::IsSparkEvent() {
    return this->isSparkEvent_;
}
int GETAnalyzer::GetNumberOfFakeEventBefore() {
    return this->nFakeEvent_;
}
int GETAnalyzer::GetNumberOfSparkEventBefore() {
    return this->nSparkEvent_;
}
bool GETAnalyzer::IsFPNChan(int chanId) {
    return (chanId == FPN_chanId_[0] || chanId == FPN_chanId_[1] || chanId == FPN_chanId_[2] || chanId == FPN_chanId_[3]) ? true : false;
}
bool GETAnalyzer::IsEvenChan(int agetId, int chanId) {
    if (this->IsFPNChan(chanId) || chanId < 0 || chanId > 67 || agetId < 0 || agetId > 3) {
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
bool GETAnalyzer::IsDeadChan(int agetId, int chanId) {
    if (chanId < 0 || chanId > 67 || agetId < 0 || agetId > 3) {
        std::cerr << "! Error in GETAnalyzer::IsDeadChan" << std::endl;
        std::cerr << "Out of range of agetId or chanId!" << std::endl;
        exit(-1);
    }
    if (this->IsFPNChan(chanId)) return false;
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
        if (this->IsEvenChan(agetId, chanId)) {  // AGET 3 even channels : only chanId 39 is dead in even channel.
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