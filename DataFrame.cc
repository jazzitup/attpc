#include "DataFrame.h"

DataFrame::DataFrame() {
    Clear();
}
DataFrame::~DataFrame() {
    PrintHeader();
}
void DataFrame::ReadHeader(std::ifstream &inFile) {
    inFile.read((char *)&(this->header_), sizeof(this->header_));
}
void DataFrame::ReadItem(std::ifstream &inFile) {
    int nItems = GetItems();
    for (int itemIdx = 0; itemIdx < nItems; itemIdx++) {
        inFile.read((char *)&(this->item_), sizeof(this->item_));
        int agetIdx = GetAgetIdx();
        int chanIdx = GetChanIdx();
        int buckIdx = GetBuckIdx();
        int sample = GetSample();
        this->ADC_[agetIdx][chanIdx][buckIdx] = sample;
    }
}
int DataFrame::GetAgetIdx() {
    return (this->item_ >> 6) & 0x03;
}
int DataFrame::GetChanIdx() {
    int chanIdx, chanIdxLowerBits, chanIdxUpperBits;
    chanIdxLowerBits = (this->item_ >> 15) & 0x01;
    chanIdxUpperBits = this->item_ & 0x3f;
    chanIdx = (chanIdxUpperBits << 1) | chanIdxLowerBits;
    return chanIdx;
}
int DataFrame::GetBuckIdx() {
    int buckIdxLowerBits, buckIdxUpperBits;
    buckIdxLowerBits = (this->item_ >> 22) & 0x03;
    buckIdxUpperBits = (this->item_ >> 8) & 0x7f;
    return (buckIdxUpperBits << 2) | buckIdxLowerBits;
}
int DataFrame::GetSample() {
    int sampleLowerBits, sampleUpperBits;
    sampleLowerBits = (this->item_ >> 24) & 0xff;
    sampleUpperBits = (this->item_ >> 16) & 0x0f;
    return (sampleUpperBits << 8) | sampleLowerBits;
}
int DataFrame::GetEventIdx() {
    return (int)((this->header_.eventIdx[0] << 24) | (this->header_.eventIdx[1] << 16) |
                 (this->header_.eventIdx[2] << 8) | (this->header_.eventIdx[3]));
}
int DataFrame::GetItems() {
    return (int)((this->header_.nItems[0] << 24) | (this->header_.nItems[1] << 16) |
                 (this->header_.nItems[2] << 8) | this->header_.nItems[3]);
}
bool DataFrame::IsFPNChannel(int chanIdx) {
    return (chanIdx == 11 || chanIdx == 22 || chanIdx == 45 || chanIdx == 56);
}
int DataFrame::GetADC(int agetIdx, int chanIdx, int buckIdx) {
    return this->ADC_[agetIdx][chanIdx][buckIdx];
}
void DataFrame::Clear() {
    for (int agetIdx = 0; agetIdx < 4; agetIdx++) {
        for (int chanIdx = 0; chanIdx < 68; chanIdx++) {
            for (int buckIdx = 0; buckIdx < 512; buckIdx++) {
                this->ADC_[agetIdx][chanIdx][buckIdx] = 0;
            }
        }
    }
}
void DataFrame::PrintHeader() {
    std::cout << std::endl;
    std::cout << "***************Header Status********************" << std::endl;
    std::cout << "MetaType : " << pow(2, (int)this->header_.metaType) << " bytes" << std::endl;
    std::cout << "FrameSize : " << (int)((this->header_.frameSize[0] << 16) | (this->header_.frameSize[1] << 8) | (this->header_.frameSize[2])) * 64 << " bytes" << std::endl;
    std::cout << "frameItemSource : " << (int)this->header_.frameItemSource << std::endl;
    std::cout << "FrameType : " << (int)this->header_.frameType[1] << std::endl;
    std::cout << "Revision : " << (int)this->header_.revision << std::endl;
    std::cout << "HeaderSize : " << (int)((this->header_.headerSize[0] << 8) | this->header_.headerSize[1]) * 64 << " bytes" << std::endl;
    std::cout << "ItemSize : " << (int)((this->header_.itemSize[0] << 8) | this->header_.itemSize[1]) * 8 << " bytes" << std::endl;
    std::cout << "nItems : " << (int)((this->header_.nItems[0] << 24) | (this->header_.nItems[1] << 16) | (this->header_.nItems[2] << 8) | (this->header_.nItems[3])) << std::endl;
    std::cout << "Cobo Idx : " << (int)this->header_.coboIdx << std::endl;
    std::cout << "AsAd Idx : " << (int)this->header_.asadIdx << std::endl;
    std::cout << "ReadOffset : " << (int)((this->header_.readOffset[0] << 8) | this->header_.readOffset[1]) << std::endl;
    std::cout << "Status : " << (int)this->header_.status << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << std::endl;
}