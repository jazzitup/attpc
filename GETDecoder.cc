#include "GETDecoder.h"

GETDecoder::GETDecoder() {
}
void GETDecoder::Open(const std::string &fileName) {
    this->grawFile_.open(fileName.c_str(), std::ios::binary);
    if (!this->grawFile_.is_open()) {
        std::cerr << "! Error in GETDecoder::Open" << std::endl;
        std::cerr << "File " << fileName.c_str() << " does not exist!" << std::endl;
        exit(-1);
    }
    this->year_ = std::stoi(fileName.substr(fileName.length() - 34, 4));
    this->month_ = std::stoi(fileName.substr(fileName.length() - 29, 2));
    this->hour_ = std::stoi(fileName.substr(fileName.length() - 26, 2));
    this->day_ = std::stoi(fileName.substr(fileName.length() - 23, 2));
    this->minute_ = std::stoi(fileName.substr(fileName.length() - 20, 2));
    this->second_ = std::stof(fileName.substr(fileName.length() - 17, 6));
    std::cout << "File open: " << fileName.c_str() << std::endl;
}
void GETDecoder::OpenFromList(const std::string &fileListName) {
    std::ifstream listFile_;
    listFile_.open(fileListName.c_str());
    if (!listFile_.is_open()) {
        std::cerr << "! Error in GETDecoder::OpenFromList" << std::endl;
        std::cerr << "File " << fileListName.c_str() << " does not exist!" << std::endl;
        exit(-1);
    }
    while (!listFile_.eof()) {
        std::string fileName;
        std::getline(listFile_, fileName);
        if (listFile_.eof()) break;
        this->fileQueue_.push(fileName);
    }
}
//
void GETDecoder::ReadHeader() {
    this->grawFile_.read((char *)&(this->header_), sizeof(Header));
}
void GETDecoder::ReadItem() {
    for (int agetId = 0; agetId < 4; agetId++) {
        for (int chanId = 0; chanId < 68; chanId++) {
            this->maxADC_[agetId][chanId] = 0;
            this->minADC_[agetId][chanId] = 4095;
        }
    }
    for (int itemId = 0; itemId < this->GetItems(); itemId++) {
        this->grawFile_.read((char *)&(this->item_), sizeof(uint32_t));
        int agetId = this->GetAgetId();
        int chanId = this->GetChanId();
        int buckId = this->GetBuckId();
        int sample = this->GetSample();
        this->ADC_[agetId][chanId][buckId] = sample;
        if (sample > this->maxADC_[agetId][chanId]) {
            this->maxADC_[agetId][chanId] = sample;
            this->timeToMax_[agetId][chanId] = buckId;
        }
        if (sample < this->minADC_[agetId][chanId]) {
            this->minADC_[agetId][chanId] = sample;
            this->timeToMin_[agetId][chanId] = buckId;
        }
    }
}
bool GETDecoder::Run() {
    if (!this->grawFile_.is_open() && this->fileQueue_.empty()) {
        std::cerr << "! Error in GETDecoder::Run" << std::endl;
        std::cerr << "File is not opened!" << std::endl;
        exit(-1);
    } else if (!this->grawFile_.is_open() && !this->fileQueue_.empty()) {
        this->Open(this->fileQueue_.front());
        this->fileQueue_.pop();
    }
    this->ReadHeader();
    if (this->grawFile_.eof() && this->fileQueue_.empty()) {
        this->grawFile_.close();
        std::cout << "[End Of Decoding]" << std::endl;
        return false;
    } else if (this->grawFile_.eof() && !this->fileQueue_.empty()) {
        this->grawFile_.close();
        std::cout << "[End Of File]" << std::endl;
        this->Run();
        return true;
    }
    this->ReadItem();
    if (this->GetEventId() == 0) {
        this->prevTime_ = this->GetEventTime();
        this->diffTime_ = 0;
    } else {
        this->diffTime_ = this->GetEventTime() - this->prevTime_;
        this->prevTime_ = this->GetEventTime();
    }
    return true;
}
// Get from header
int GETDecoder::GetEventId() {
    return ((this->header_.eventId[0] << 24) | (this->header_.eventId[1] << 16) |
            (this->header_.eventId[2] << 8) | (this->header_.eventId[3]));
}
uint64_t GETDecoder::GetEventTime() {
    return (((uint64_t)this->header_.eventTime[0] << 40) | ((uint64_t)this->header_.eventTime[1] << 32) |
            ((uint64_t)this->header_.eventTime[2] << 24) | ((uint64_t)this->header_.eventTime[3] << 16) |
            ((uint64_t)this->header_.eventTime[4] << 8) | ((uint64_t)this->header_.eventTime[5]));
}
uint64_t GETDecoder::GetDiffTime() {
    return this->diffTime_;
}
// Get from item
// Frame Item Format : aacc cccc | cbbb bbbb | bb00 ssss | ssss ssss
// Read in reverse order : ssss ssss | bb00 ssss | cbbb bbbb | aacc cccc
// a : agetId, c : chanId, b : buckId, s : sample
int GETDecoder::GetAgetId() {
    return (this->item_ >> 6) & 0x03;
}
int GETDecoder::GetChanId() {
    int chanId, chanIdLowerBits, chanIdUpperBits;
    chanIdLowerBits = (this->item_ >> 15) & 0x01;
    chanIdUpperBits = this->item_ & 0x3f;
    chanId = (chanIdUpperBits << 1) | chanIdLowerBits;
    return chanId;
}
int GETDecoder::GetBuckId() {
    int buckIdLowerBits, buckIdUpperBits;
    buckIdLowerBits = (this->item_ >> 22) & 0x03;
    buckIdUpperBits = (this->item_ >> 8) & 0x7f;
    return (buckIdUpperBits << 2) | buckIdLowerBits;
}
int GETDecoder::GetSample() {
    int sampleLowerBits, sampleUpperBits;
    sampleLowerBits = (this->item_ >> 24) & 0xff;
    sampleUpperBits = (this->item_ >> 16) & 0x0f;
    return (sampleUpperBits << 8) | sampleLowerBits;
}
int GETDecoder::GetItems() {
    return (int)((this->header_.nItems[0] << 24) | (this->header_.nItems[1] << 16) |
                 (this->header_.nItems[2] << 8) | this->header_.nItems[3]);
}
//
int GETDecoder::GetADC(int agetId, int chanId, int buckId) {
    return this->ADC_[agetId][chanId][buckId];
}
int GETDecoder::GetMaxADC(int agetId, int chanId) {
    return this->maxADC_[agetId][chanId];
}
int GETDecoder::GetMinADC(int agetId, int chanId) {
    return this->minADC_[agetId][chanId];
}
int GETDecoder::GetTimeToMax(int agetId, int chanId) {
    return this->timeToMax_[agetId][chanId];
}
int GETDecoder::GetTimeToMin(int agetId, int chanId) {
    return this->timeToMin_[agetId][chanId];
}